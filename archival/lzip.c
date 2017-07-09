/*
 * lzip implementation for busybox
 *
 * Copyright (C) 2012-2016 Antonio Diaz Diaz.
 *
 * Licensed under GPLv2 or later, see file LICENSE in this source tree.
 */

//config:config LZIP
//config:	bool "lzip"
//config:	default y
//config:	help
//config:	  Lzip is a lossless data compressor with a user interface similar to
//config:	  the one of gzip or bzip2. Lzip is about as fast as gzip, compresses
//config:	  most files more than bzip2, and is better than both from a data
//config:	  recovery perspective.

//applet:IF_LZIP(APPLET(lzip, BB_DIR_USR_BIN, BB_SUID_DROP))
//kbuild:lib-$(CONFIG_LZIP) += lzip.o bbunzip.o

//usage:#define lzip_trivial_usage
//usage:       "[-123456789c"
//usage:	IF_LUNZIP("d") "f"
//usage:	IF_LUNZIP("t")
//usage:       "] [-m MATCH_LENGTH] [-s DICT_SIZE] [FILE]..."
//usage:#define lzip_full_usage "\n\n"
//usage:       "Compress FILEs (or stdin) with lzip algorithm\n"
//usage:     "\n	-1..9	Compression level"
//usage:     "\n	-c	Write to stdout"
//usage:	IF_LUNZIP("\n	-d	Decompress")
//usage:     "\n	-f	Force"
//usage:     "\n	-m	Match length limit [36]"
//usage:     "\n	-s	Dictionary size limit [8MiB]"
//usage:	IF_LUNZIP("\n	-t	Test compressed file integrity")


#include "libbb.h"
#include "bb_archive.h"
#include "libarchive/lzip.h"


#if CHAR_BIT != 8
#error "Environments where CHAR_BIT != 8 are not supported."
#endif


static void CRC32_update_byte(uint32_t * crc, const uint8_t byte)
{
	*crc = global_crc32_table[(*crc ^ byte) & 0xFF] ^ (*crc >> 8);
}


enum { max_num_trials = 1 << 12,
	price_shift_bits = 6
};


static uint8_t * dis_slots;

static void Dis_slots_init(void)
{
	int i, size, slot;
	dis_slots = xmalloc((1 << 10) * sizeof dis_slots[0]);

	for (slot = 0; slot < 4; ++slot) dis_slots[slot] = slot;
	for (i = 4, size = 2, slot = 4; slot < 20; slot += 2) {
		memset(&dis_slots[i], slot, size);
		memset(&dis_slots[i + size], slot + 1, size);
		size <<= 1;
		i += size;
	}
}

static uint8_t get_slot(const unsigned dis)
{
	if (dis < (1 << 10)) return dis_slots[dis];
	if (dis < (1 << 19)) return dis_slots[dis>> 9] + 18;
	if (dis < (1 << 28)) return dis_slots[dis>>18] + 36;
	return dis_slots[dis>>27] + 54;
}


static int * prob_prices;

static void Prob_prices_init(void)
{
	const int num_bits = (bit_model_total_bits - 2);
	int i, j = 1, end = 2;
	prob_prices = xmalloc((bit_model_total >> 2) * sizeof prob_prices[0]);

	prob_prices[0] = bit_model_total_bits << price_shift_bits;
	for (i = num_bits - 1; i >= 0; --i, end <<= 1) {
		for (; j < end; ++j)
			prob_prices[j] = (i << price_shift_bits) +
				(((end - j) << price_shift_bits) >> (num_bits - i - 1));
	}
}

static inline int get_price(const int probability)
{
	return prob_prices[probability >> 2];
}


static inline int price0(const Bit_model probability)
{
	return get_price(probability);
}

static inline int price1(const Bit_model probability)
{
	return get_price(bit_model_total - probability);
}

static int price_bit(const Bit_model bm, const int bit)
{
	if (bit) return price1(bm);
	else return price0(bm);
}


static int price_symbol(const Bit_model bm[], int symbol,
						const int num_bits)
{
	int price = 0;
	symbol |= (1 << num_bits);
	while (symbol > 1) {
		const int bit = symbol & 1;
		symbol >>= 1;
		price += price_bit(bm[symbol], bit);
	}
	return price;
}


static int price_symbol_reversed(const Bit_model bm[], int symbol,
					const int num_bits)
{
	int price = 0;
	int model = 1;
	int i;
	for (i = num_bits; i > 0; --i) {
		const int bit = symbol & 1;
		price += price_bit(bm[model], bit);
		model = (model << 1) | bit;
		symbol >>= 1;
	}
	return price;
}


static int price_matched(const Bit_model bm[], int symbol, int match_byte)
{
	int price = 0;
	int mask = 0x100;
	symbol |= mask;

	do {
		int match_bit, bit;
		match_byte <<= 1;
		match_bit = match_byte & mask;
		symbol <<= 1;
		bit = symbol & 0x100;
		price += price_bit( bm[match_bit+(symbol>>9)+mask], bit );
		mask &= ~(match_byte ^ symbol);	/* if( match_bit != bit ) mask = 0; */
	}
	while( symbol < 0x10000 );
	return price;
}


enum {	/* bytes to keep in buffer before dictionary */
	before_size = max_num_trials + 1,
	/* bytes to keep in buffer after pos */
	after_size = max_match_len,
	num_prev_positions4 = 1 << 20,
	num_prev_positions3 = 1 << 18,
	num_prev_positions2 = 1 << 16,
	num_prev_positions = num_prev_positions4 + num_prev_positions3 +
		num_prev_positions2
};

struct Matchfinder {
	unsigned long long partial_data_pos;
	uint8_t *buffer;		/* input buffer */
	int32_t *prev_positions;	/* last seen position of key */
	int32_t *prev_pos_tree;		/* previous positions of key */
	int match_len_limit;
	int buffer_size;
	int dictionary_size;	/* bytes to keep in buffer before pos */
	int pos;		/* current pos in buffer */
	int cyclic_pos;		/* current pos in dictionary */
	int stream_pos;		/* first byte not yet read from file */
	int pos_limit;		/* when reached, a new block must be read */
	int cycles;
	bool at_stream_end;	/* stream_pos shows real end of file */
};

static bool Mf_read_block(struct Matchfinder *const mf)
{
	if (!mf->at_stream_end && mf->stream_pos < mf->buffer_size) {
		const int size = mf->buffer_size - mf->stream_pos;
		const int rd = full_read(STDIN_FILENO,
					mf->buffer + mf->stream_pos, size);
		mf->stream_pos += rd;
		if (rd < size) {
			mf->at_stream_end = true;
			mf->pos_limit = mf->buffer_size;
		}
	}
	return mf->pos < mf->stream_pos;
}

static void Mf_normalize_pos(struct Matchfinder *const mf)
{
	if (!mf->at_stream_end) {
		int i;
		const int offset = mf->pos - mf->dictionary_size - before_size;
		const int size = mf->stream_pos - offset;
		memmove(mf->buffer, mf->buffer + offset, size);
		mf->partial_data_pos += offset;
		mf->pos -= offset;
		mf->stream_pos -= offset;
		for (i = 0; i < num_prev_positions; ++i)
			if (mf->prev_positions[i] >= 0)
				mf->prev_positions[i] -= offset;
		for (i = 0; i < 2 * mf->dictionary_size; ++i)
			if (mf->prev_pos_tree[i] >= 0)
				mf->prev_pos_tree[i] -= offset;
		Mf_read_block(mf);
	}
}

static bool Mf_init(struct Matchfinder *const mf, const int dict_size,
			const int match_len_limit)
{
	const int buffer_size_limit = (2 * dict_size) + before_size + after_size;
	int i;

	mf->partial_data_pos = 0;
	mf->match_len_limit = match_len_limit;
	mf->prev_positions = (int32_t *) malloc(num_prev_positions * sizeof(int32_t));
	if (!mf->prev_positions) return false;
	mf->pos = 0;
	mf->cyclic_pos = 0;
	mf->stream_pos = 0;
	mf->cycles = (match_len_limit < max_match_len) ?
		16 + (match_len_limit / 2) : 256;
	mf->at_stream_end = false;

	for (i = 0; i < num_prev_positions; ++i)
		mf->prev_positions[i] = -1;
	mf->buffer_size = MAX(65536, dict_size);
	mf->buffer = (uint8_t *) malloc(mf->buffer_size);
	if (!mf->buffer) {
		free(mf->prev_positions);
		return false;
	}
	if (Mf_read_block(mf) && !mf->at_stream_end &&
		mf->buffer_size < buffer_size_limit) {
		uint8_t *tmp;
		mf->buffer_size = buffer_size_limit;
		tmp = (uint8_t *) realloc(mf->buffer, mf->buffer_size);
		if (!tmp) {
			free(mf->buffer);
			free(mf->prev_positions);
			return false;
		}
		mf->buffer = tmp;
		Mf_read_block(mf);
	}
	if (mf->at_stream_end && mf->stream_pos < dict_size)
		mf->dictionary_size = MAX(min_dictionary_size, mf->stream_pos);
	else
		mf->dictionary_size = dict_size;
	mf->pos_limit = mf->buffer_size;
	if (!mf->at_stream_end) mf->pos_limit -= after_size;
	mf->prev_pos_tree =
		(int32_t *) malloc(2 * mf->dictionary_size * sizeof(int32_t));
	if (!mf->prev_pos_tree) {
		free(mf->buffer);
		free(mf->prev_positions);
		return false;
	}
	return true;
}

static void Mf_free(struct Matchfinder *const mf)
{
	free(mf->prev_pos_tree);
	free(mf->buffer);
	free(mf->prev_positions);
}

static inline uint8_t Mf_peek(const struct Matchfinder *const mf,
				const int distance)
{
	return mf->buffer[mf->pos-distance];
}

static inline int Mf_available_bytes(const struct Matchfinder *const mf)
{
	return mf->stream_pos - mf->pos;
}

static inline unsigned long long
Mf_data_position(const struct Matchfinder *const mf)
{
	return mf->partial_data_pos + mf->pos;
}

static inline bool Mf_finished(const struct Matchfinder *const mf)
{
	return mf->at_stream_end && mf->pos >= mf->stream_pos;
}

static inline const uint8_t *
Mf_ptr_to_current_pos(const struct Matchfinder *const mf)
{
	return mf->buffer + mf->pos;
}

static int Mf_true_match_len(const struct Matchfinder *const mf,
			const int index, const int distance, int len_limit)
{
	const uint8_t *const data = mf->buffer + mf->pos + index;
	int i = 0;
	if (index + len_limit > Mf_available_bytes(mf))
		len_limit = Mf_available_bytes(mf) - index;
	while (i < len_limit && data[i - distance] == data[i]) ++i;
	return i;
}

static void Mf_move_pos(struct Matchfinder *const mf)
{
	if (++mf->cyclic_pos >= mf->dictionary_size) mf->cyclic_pos = 0;
	if (++mf->pos >= mf->pos_limit) Mf_normalize_pos(mf);
}

static int Mf_longest_match_len(struct Matchfinder *const mf,
				int *const distances)
{
	int32_t *ptr0 = mf->prev_pos_tree + (mf->cyclic_pos << 1);
	int32_t *ptr1 = ptr0 + 1;
	int32_t *newptr;
	const uint8_t *newdata;
	int len = 0, len0 = 0, len1 = 0;
	int maxlen = min_match_len - 1;
	const int min_pos = (mf->pos >= mf->dictionary_size) ?
		(mf->pos - mf->dictionary_size + 1) : 0;
	const uint8_t *const data = mf->buffer + mf->pos;
	int count, delta, key2, key3, key4, newpos, tmp;
	int len_limit = mf->match_len_limit;

	if (len_limit > Mf_available_bytes(mf)) {
		len_limit = Mf_available_bytes(mf);
		if (len_limit < 4) return 0;
	}

	key2 = num_prev_positions4 + num_prev_positions3 +
		(((int) data[0] << 8) | data[1]);
	tmp = global_crc32_table[data[0]] ^ data[1] ^ ((uint32_t) data[2] << 8);
	key3 = num_prev_positions4 + (int) (tmp & (num_prev_positions3 - 1));
	key4 = (int) ((tmp ^ (global_crc32_table[data[3]] << 5)) &
				(num_prev_positions4 - 1));

	if (distances) {
		int np = mf->prev_positions[key2];
		if (np >= min_pos) {
			distances[2] = mf->pos - np - 1;
			maxlen = 2;
		} else
			distances[2] = 0x7FFFFFFF;
		np = mf->prev_positions[key3];
		if (np >= min_pos && mf->buffer[np] == data[0]) {
			distances[3] = mf->pos - np - 1;
			maxlen = 3;
		} else
			distances[3] = 0x7FFFFFFF;
		distances[4] = 0x7FFFFFFF;
	}

	mf->prev_positions[key2] = mf->pos;
	mf->prev_positions[key3] = mf->pos;
	newpos = mf->prev_positions[key4];
	mf->prev_positions[key4] = mf->pos;

	for (count = mf->cycles;;) {
		if (newpos < min_pos || --count < 0) {
			*ptr0 = *ptr1 = -1;
			break;
		}
		newdata = mf->buffer + newpos;
		while (len < len_limit && newdata[len] == data[len]) ++len;

		delta = mf->pos - newpos;
		if (distances)
			while (maxlen < len)
				distances[++maxlen] = delta - 1;

		newptr = mf->prev_pos_tree +
			((mf->cyclic_pos - delta +
			((mf->cyclic_pos >= delta) ? 0 : mf->dictionary_size)) << 1);

		if (len < len_limit) {
			if (newdata[len] < data[len]) {
				*ptr0 = newpos;
				ptr0 = newptr + 1;
				newpos = *ptr0;
				len0 = len;
				if (len1 < len) len = len1;
			} else {
				*ptr1 = newpos;
				ptr1 = newptr;
				newpos = *ptr1;
				len1 = len;
				if (len0 < len) len = len0;
			}
		} else {
			*ptr0 = newptr[0];
			*ptr1 = newptr[1];
			break;
		}
	}
	if (distances) {
		if (distances[3] > distances[4])
			distances[3] = distances[4];
		if (distances[2] > distances[3])
			distances[2] = distances[3];
	}
	return maxlen;
}


enum { re_buffer_size = 16384 };

struct Range_encoder {
	uint64_t low;
	unsigned long long partial_member_pos;
	uint8_t *buffer;		/* output buffer */
	int pos;			/* current pos in buffer */
	uint32_t range;
	unsigned ff_count;
	uint8_t cache;
	bool write_error;
};

static void Re_flush_data(struct Range_encoder *const renc)
{
	if (renc->pos > 0) {
		if (full_write(STDOUT_FILENO, renc->buffer, renc->pos) != renc->pos)
			renc->write_error = true;
		renc->partial_member_pos += renc->pos;
		renc->pos = 0;
	}
}

static void Re_put_byte(struct Range_encoder *const renc, const uint8_t b)
{
	renc->buffer[renc->pos] = b;
	if (++renc->pos >= re_buffer_size) Re_flush_data(renc);
}

static void Re_shift_low(struct Range_encoder *const renc)
{
	const bool carry = (renc->low > 0xFFFFFFFFU);
	if (carry || renc->low < 0xFF000000U) {
		Re_put_byte(renc, renc->cache + carry);
		for (; renc->ff_count > 0; --renc->ff_count)
			Re_put_byte(renc, 0xFF + carry);
		renc->cache = renc->low >> 24;
	} else
		++renc->ff_count;
	renc->low = (renc->low & 0x00FFFFFFU) << 8;
}

static bool Re_init(struct Range_encoder *const renc)
{
	renc->low = 0;
	renc->partial_member_pos = 0;
	renc->buffer = (uint8_t *) malloc(re_buffer_size);
	if (!renc->buffer) return false;
	renc->pos = 0;
	renc->range = 0xFFFFFFFFU;
	renc->ff_count = 0;
	renc->cache = 0;
	renc->write_error = false;
	return true;
}

static inline void Re_free(struct Range_encoder *const renc)
{
	free(renc->buffer);
}

static inline unsigned long long
Re_member_position(const struct Range_encoder *const renc)
{
	return renc->partial_member_pos + renc->pos + renc->ff_count;
}

static void Re_flush(struct Range_encoder *const renc)
{
	int i;
	for (i = 0; i < 5; ++i) Re_shift_low(renc);
}

static void Re_encode(struct Range_encoder *const renc,
				const int symbol, const int num_bits)
{
	int i;
	for (i = num_bits - 1; i >= 0; --i) {
		renc->range >>= 1;
		if ((symbol >> i) & 1) renc->low += renc->range;
		if (renc->range <= 0x00FFFFFFU) {
			renc->range <<= 8;
			Re_shift_low(renc);
		}
	}
}

static void Re_encode_bit(struct Range_encoder *const renc,
				Bit_model * const probability, const int bit)
{
	const uint32_t bound = (renc->range >> bit_model_total_bits) * *probability;
	if (!bit) {
		renc->range = bound;
		*probability += (bit_model_total - *probability) >> bit_model_move_bits;
	} else {
		renc->low += bound;
		renc->range -= bound;
		*probability -= *probability >> bit_model_move_bits;
	}
	if (renc->range <= 0x00FFFFFFU) {
		renc->range <<= 8;
		Re_shift_low(renc);
	}
}

static void Re_encode_tree(struct Range_encoder *const renc,
				Bit_model bm[], const int symbol,
				const int num_bits)
{
	int mask = (1 << (num_bits - 1));
	int model = 1;
	int i;
	for (i = num_bits; i > 0; --i, mask >>= 1) {
		const int bit = (symbol & mask);
		Re_encode_bit(renc, &bm[model], bit);
		model <<= 1;
		if (bit) model |= 1;
	}
}

static void Re_encode_tree_reversed(struct Range_encoder *const renc,
			Bit_model bm[], int symbol, const int num_bits)
{
	int model = 1;
	int i;
	for (i = num_bits; i > 0; --i) {
		const int bit = symbol & 1;
		Re_encode_bit(renc, &bm[model], bit);
		model = (model << 1) | bit;
		symbol >>= 1;
	}
}

static void Re_encode_matched(struct Range_encoder *const renc,
				Bit_model bm[], int symbol, int match_byte)
{
	int mask = 0x100;
	symbol |= mask;

	do {
		int match_bit, bit;
		match_byte <<= 1;
		match_bit = match_byte & mask;
		symbol <<= 1;
		bit = symbol & 0x100;
		Re_encode_bit( renc, &bm[match_bit+(symbol>>9)+mask], bit );
		mask &= ~(match_byte ^ symbol);	/* if( match_bit != bit ) mask = 0; */
	}
	while( symbol < 0x10000 );
}

static void Re_encode_len( struct Range_encoder * const renc,
				struct Len_model * const lm,
				int symbol, const int pos_state )
{
	bool bit = ( ( symbol -= min_match_len ) >= len_low_symbols );
	Re_encode_bit( renc, &lm->choice1, bit );
	if( !bit )
		Re_encode_tree( renc, lm->bm_low[pos_state], symbol, len_low_bits );
	else {
		bit = ( symbol >= len_low_symbols + len_mid_symbols );
		Re_encode_bit( renc, &lm->choice2, bit );
		if( !bit )
			Re_encode_tree( renc, lm->bm_mid[pos_state],
				symbol - len_low_symbols, len_mid_bits );
		else
			Re_encode_tree( renc, lm->bm_high,
			symbol - len_low_symbols - len_mid_symbols, len_high_bits );
	}
}


struct Len_encoder {
	struct Len_model lm;
	int len_symbols;
	int prices[pos_states][max_len_symbols];
	int counters[pos_states];
};

static void Lee_update_prices(struct Len_encoder *const le, const int pos_state)
{
	int *const pps = le->prices[pos_state];
	int tmp = price0(le->lm.choice1);
	int len = 0;

	for (; len < len_low_symbols && len < le->len_symbols; ++len)
		pps[len] = tmp +
			price_symbol(le->lm.bm_low[pos_state], len, len_low_bits);
	tmp = price1(le->lm.choice1);
	for (; len < len_low_symbols + len_mid_symbols && len < le->len_symbols; ++len)
		pps[len] = tmp + price0(le->lm.choice2) +
			price_symbol(le->lm.bm_mid[pos_state],
					len - len_low_symbols, len_mid_bits);
	for (; len < le->len_symbols; ++len)
		/* using 4 slots per value makes "Lee_price" faster */
		le->prices[3][len] = le->prices[2][len] =
		le->prices[1][len] = le->prices[0][len] =
			tmp + price1(le->lm.choice2) +
			price_symbol(le->lm.bm_high,
				len - len_low_symbols - len_mid_symbols,
				len_high_bits);
	le->counters[pos_state] = le->len_symbols;
}

static void Lee_init(struct Len_encoder *const le, const int len_limit)
{
	int i;
	Lm_init(&le->lm);
	le->len_symbols = len_limit + 1 - min_match_len;
	for (i = 0; i < pos_states; ++i) Lee_update_prices(le, i);
}

static void Lee_encode(struct Len_encoder *const le,
			struct Range_encoder *const renc,
			int symbol, const int pos_state)
{
	Re_encode_len(renc, &le->lm, symbol, pos_state);
	if (--le->counters[pos_state] <= 0)
		Lee_update_prices(le, pos_state);
}

static int Lee_price(const struct Len_encoder *const le,
			const int symbol, const int pos_state)
{
	return le->prices[pos_state][symbol - min_match_len];
}


enum { infinite_price = 0x0FFFFFFF,
	num_rep_distances = 4			/* must be 4 */
};

struct Trial {
	State state;
	int price;	/* dual use var; cumulative price, match length */
	int dis;	/* rep index or match distance. (-1 for literal) */
	int prev_index;	/* index of prev trial in trials[] */
	int reps[num_rep_distances];
};

static void Tr_update(struct Trial *const trial, const int pr,
			const int distance, const int p_i)
{
	if (pr < trial->price) {
		trial->price = pr;
		trial->dis = distance;
		trial->prev_index = p_i;
	}
}


struct LZ_encoder {
	int longest_match_found;
	uint32_t crc;

	Bit_model bm_literal[1<<literal_context_bits][0x300];
	Bit_model bm_match[states][pos_states];
	Bit_model bm_rep[states];
	Bit_model bm_rep0[states];
	Bit_model bm_rep1[states];
	Bit_model bm_rep2[states];
	Bit_model bm_len[states][pos_states];
	Bit_model bm_dis_slot[len_states][1<<dis_slot_bits];
	Bit_model bm_dis[modeled_distances-end_dis_model];
	Bit_model bm_align[dis_align_size];

	struct Matchfinder *matchfinder;
	struct Range_encoder renc;
	struct Len_encoder match_len_encoder;
	struct Len_encoder rep_len_encoder;

	int match_distances[max_match_len+1];
	struct Trial trials[max_num_trials];

	int dis_slot_prices[len_states][2*max_dictionary_bits];
	int dis_prices[len_states][modeled_distances];
	int align_prices[dis_align_size];
	int align_price_count;
	int num_dis_slots;
};

static void LZe_fill_align_prices(struct LZ_encoder *const e)
{
	int i;
	for (i = 0; i < dis_align_size; ++i)
		e->align_prices[i] =
			price_symbol_reversed(e->bm_align, i, dis_align_bits);
	e->align_price_count = dis_align_size;
}

static bool LZe_init(struct LZ_encoder *const e,
			struct Matchfinder *const mf, const File_header header)
{
	int i;
	e->longest_match_found = 0;
	e->crc = 0xFFFFFFFFU;
	Bm_array_init(&e->bm_literal[0][0], lz_num_models);
	e->matchfinder = mf;
	if (!Re_init(&e->renc)) return false;
	Lee_init(&e->match_len_encoder, e->matchfinder->match_len_limit);
	Lee_init(&e->rep_len_encoder, e->matchfinder->match_len_limit);
	LZe_fill_align_prices(e);
	e->num_dis_slots = 2 * real_bits(e->matchfinder->dictionary_size - 1);
	for (i = 0; i < Fh_size; ++i)
		Re_put_byte(&e->renc, header[i]);
	return true;
}

static inline void LZe_free(struct LZ_encoder *const e)
{
	Re_free(&e->renc);
}

static inline unsigned LZe_crc(const struct LZ_encoder *const e)
{
	return e->crc ^ 0xFFFFFFFFU;
}

	/* move-to-front dis in/into reps if( dis > 0 ) */
static void mtf_reps(const int dis, int reps[num_rep_distances])
{
	int i;
	if (dis >= num_rep_distances) {
		for (i = num_rep_distances - 1; i > 0; --i)
			reps[i] = reps[i - 1];
		reps[0] = dis - num_rep_distances;
	} else if (dis > 0) {
		const int distance = reps[dis];
		for (i = dis; i > 0; --i)
			reps[i] = reps[i - 1];
		reps[0] = distance;
	}
}

static int LZe_price_shortrep(const struct LZ_encoder *const e,
				const State state, const int pos_state)
{
	return price0(e->bm_rep0[state]) + price0(e->bm_len[state][pos_state]);
}

static int LZe_price_rep(const struct LZ_encoder *const e, const int rep,
				const State state, const int pos_state)
{
	int price;
	if (rep == 0)
		return price0(e->bm_rep0[state]) +
			price1(e->bm_len[state][pos_state]);
	price = price1(e->bm_rep0[state]);
	if (rep == 1)
		price += price0(e->bm_rep1[state]);
	else {
		price += price1(e->bm_rep1[state]);
		price += price_bit(e->bm_rep2[state], rep - 2);
	}
	return price;
}

static int LZe_price_dis(const struct LZ_encoder *const e,
				const int dis, const int len_state)
{
	if (dis < modeled_distances)
		return e->dis_prices[len_state][dis];
	else
		return e->dis_slot_prices[len_state][get_slot(dis)] +
			e->align_prices[dis & (dis_align_size - 1)];
}

static int LZe_price_pair(const struct LZ_encoder *const e,
				const int dis, const int len,
				const int pos_state)
{
	if (len <= min_match_len && dis >= modeled_distances)
		return infinite_price;
	return Lee_price(&e->match_len_encoder, len, pos_state) +
		LZe_price_dis(e, dis, get_len_state(len));
}

static int LZe_price_literal(const struct LZ_encoder *const e,
				uint8_t prev_byte, uint8_t symbol)
{
	return price_symbol(e->bm_literal[get_lit_state(prev_byte)], symbol, 8);
}

static int LZe_price_matched(const struct LZ_encoder *const e,
				uint8_t prev_byte, uint8_t symbol,
				uint8_t match_byte)
{
	return price_matched(e->bm_literal[get_lit_state(prev_byte)], symbol,
				match_byte);
}

static void LZe_encode_literal(struct LZ_encoder *const e,
				uint8_t prev_byte, uint8_t symbol)
{
	Re_encode_tree(&e->renc,
		e->bm_literal[get_lit_state(prev_byte)], symbol, 8);
}

static void LZe_encode_matched(struct LZ_encoder *const e,
				uint8_t prev_byte, uint8_t symbol,
				uint8_t match_byte)
{
	Re_encode_matched(&e->renc, e->bm_literal[get_lit_state(prev_byte)],
			symbol, match_byte);
}

static void LZe_encode_pair(struct LZ_encoder *const e,
				const unsigned dis, const int len,
				const int pos_state)
{
	const int dis_slot = get_slot(dis);
	Lee_encode(&e->match_len_encoder, &e->renc, len, pos_state);
	Re_encode_tree(&e->renc, e->bm_dis_slot[get_len_state(len)], dis_slot,
			dis_slot_bits);

	if (dis_slot >= start_dis_model) {
		const int direct_bits = (dis_slot >> 1) - 1;
		const unsigned base = (2 | (dis_slot & 1)) << direct_bits;
		const unsigned direct_dis = dis - base;

		if (dis_slot < end_dis_model)
			Re_encode_tree_reversed(&e->renc,
					e->bm_dis + base - dis_slot - 1,
					direct_dis, direct_bits);
		else {
			Re_encode(&e->renc, direct_dis >> dis_align_bits,
					direct_bits - dis_align_bits);
			Re_encode_tree_reversed(&e->renc, e->bm_align,
					direct_dis, dis_align_bits);
			if (--e->align_price_count <= 0)
				LZe_fill_align_prices(e);
		}
	}
}

static int LZe_read_match_distances(struct LZ_encoder *const e)
{
	int len = Mf_longest_match_len(e->matchfinder, e->match_distances);
	if (len == e->matchfinder->match_len_limit && len < max_match_len)
		len += Mf_true_match_len(e->matchfinder, len,
					e->match_distances[len] + 1,
					max_match_len - len);
	return len;
}

static void LZe_move_pos(struct LZ_encoder *const e, int n)
{
	while (true) {
		Mf_move_pos(e->matchfinder);
		if( --n <= 0 ) break;
		Mf_longest_match_len(e->matchfinder, 0);
	}
}

static void LZe_backward(struct LZ_encoder *const e, int cur)
{
	int *const dis = &e->trials[cur].dis;
	while (cur > 0) {
		const int prev_index = e->trials[cur].prev_index;
		struct Trial *const prev_trial = &e->trials[prev_index];
		prev_trial->price = cur - prev_index;	/* len */
		cur = *dis;
		*dis = prev_trial->dis;
		prev_trial->dis = cur;
		cur = prev_index;
	}
}

	/* End Of Stream mark => (dis == 0xFFFFFFFFU, len == min_match_len) */
static void LZe_full_flush(struct LZ_encoder *const e, const State state)
{
	int i;
	const int pos_state = Mf_data_position(e->matchfinder) & pos_state_mask;
	File_trailer trailer;
	Re_encode_bit(&e->renc, &e->bm_match[state][pos_state], 1);
	Re_encode_bit(&e->renc, &e->bm_rep[state], 0);
	LZe_encode_pair(e, 0xFFFFFFFFU, min_match_len, pos_state);
	Re_flush(&e->renc);
	Ft_set_data_crc(trailer, LZe_crc(e));
	Ft_set_data_size(trailer, Mf_data_position(e->matchfinder));
	Ft_set_member_size(trailer, Re_member_position(&e->renc) + Ft_size);
	for (i = 0; i < Ft_size; ++i)
		Re_put_byte(&e->renc, trailer[i]);
	Re_flush_data(&e->renc);
}


static void LZe_update_distance_prices(struct LZ_encoder *const e)
{
	int dis, len_state;
	for (dis = start_dis_model; dis < modeled_distances; ++dis) {
		const int dis_slot = dis_slots[dis];
		const int direct_bits = (dis_slot >> 1) - 1;
		const int base = (2 | (dis_slot & 1)) << direct_bits;
		const int price =
			price_symbol_reversed(e->bm_dis + base - dis_slot - 1,
						dis - base, direct_bits);
		for (len_state = 0; len_state < len_states; ++len_state)
			e->dis_prices[len_state][dis] = price;
	}

	for (len_state = 0; len_state < len_states; ++len_state) {
		int *const dsp = e->dis_slot_prices[len_state];
		int *const dp = e->dis_prices[len_state];
		const Bit_model *const bmds = e->bm_dis_slot[len_state];
		int slot = 0;
		for (; slot < end_dis_model; ++slot)
			dsp[slot] = price_symbol(bmds, slot, dis_slot_bits);
		for (; slot < e->num_dis_slots; ++slot)
			dsp[slot] = price_symbol(bmds, slot, dis_slot_bits) +
				((((slot >> 1) - 1) - dis_align_bits) << price_shift_bits);

		for (dis = 0; dis < start_dis_model; ++dis)
			dp[dis] = dsp[dis];
		for (; dis < modeled_distances; ++dis)
			dp[dis] += dsp[dis_slots[dis]];
	}
}


/* Returns the number of bytes advanced (ahead).
   trials[0]..trials[ahead-1] contain the steps to encode.
   ( trials[0].dis == -1 && trials[0].price == 1 ) means literal.
   A match/rep longer or equal than match_len_limit finishes the sequence.
*/
static int LZe_sequence_optimizer(struct LZ_encoder *const e,
				const int reps[num_rep_distances],
				const State state)
{
	int main_len, i, rep, cur = 0, num_trials;
	int replens[num_rep_distances];
	int rep_index = 0;

	if (e->longest_match_found > 0) {	/* from previous call */
		main_len = e->longest_match_found;
		e->longest_match_found = 0;
	} else
		main_len = LZe_read_match_distances(e);

	for (i = 0; i < num_rep_distances; ++i) {
		replens[i] = Mf_true_match_len(e->matchfinder, 0, reps[i] + 1,
						max_match_len);
		if (replens[i] > replens[rep_index]) rep_index = i;
	}
	if (replens[rep_index] >= e->matchfinder->match_len_limit) {
		e->trials[0].dis = rep_index;
		e->trials[0].price = replens[rep_index];
		LZe_move_pos(e, replens[rep_index]);
		return replens[rep_index];
	}

	if (main_len >= e->matchfinder->match_len_limit) {
		e->trials[0].dis =
			e->match_distances[e->matchfinder->match_len_limit] +
			num_rep_distances;
		e->trials[0].price = main_len;
		LZe_move_pos(e, main_len);
		return main_len;
	}

	{
	const int pos_state = Mf_data_position(e->matchfinder) & pos_state_mask;
	const int match_price = price1(e->bm_match[state][pos_state]);
	const int rep_match_price = match_price + price1(e->bm_rep[state]);
	const uint8_t prev_byte = Mf_peek(e->matchfinder, 1);
	const uint8_t cur_byte = Mf_peek(e->matchfinder, 0);
	const uint8_t match_byte = Mf_peek(e->matchfinder, reps[0] + 1);

	e->trials[0].state = state;
	for (i = 0; i < num_rep_distances; ++i)
		e->trials[0].reps[i] = reps[i];
	e->trials[1].dis = -1;				/* literal */
	e->trials[1].prev_index = 0;
	e->trials[1].price = price0(e->bm_match[state][pos_state]);
	if (St_is_char(state))
		e->trials[1].price +=
			LZe_price_literal(e, prev_byte, cur_byte);
	else
		e->trials[1].price +=
			LZe_price_matched(e, prev_byte, cur_byte, match_byte);

	if (match_byte == cur_byte)
		Tr_update(&e->trials[1], rep_match_price +
			LZe_price_shortrep(e, state, pos_state), 0, 0);

	if (main_len < min_match_len) {
		e->trials[0].dis = e->trials[1].dis;
		e->trials[0].price = 1;
		Mf_move_pos(e->matchfinder);
		return 1;
	}

	if (main_len <= replens[rep_index]) {
		int len;

		main_len = replens[rep_index];
		for (len = min_match_len; len <= main_len; ++len)
			e->trials[len].price = infinite_price;
	} else {
		int len;
		const int normal_match_price =
			match_price + price0(e->bm_rep[state]);
		for (len = min_match_len; len <= main_len; ++len) {
			e->trials[len].dis =
				e->match_distances[len] + num_rep_distances;
			e->trials[len].prev_index = 0;
			e->trials[len].price = normal_match_price +
				LZe_price_pair(e, e->match_distances[len],
							len, pos_state);
		}
	}

	for (rep = 0; rep < num_rep_distances; ++rep) {
		const int price = rep_match_price +
			LZe_price_rep(e, rep, state, pos_state);
		int len;
		for (len = min_match_len; len <= replens[rep]; ++len)
			Tr_update(&e->trials[len], price +
				Lee_price(&e->rep_len_encoder, len, pos_state),
				rep, 0);
	}
	}

	num_trials = main_len;

	while (true) {				/* price optimization loop */
		struct Trial *cur_trial, *next_trial;
		int newlen, pos_state, prev_index, len_limit;
		int next_price, match_price, rep_match_price;
		uint8_t prev_byte, cur_byte, match_byte;

		Mf_move_pos(e->matchfinder);
		if (++cur >= num_trials) {	/* no more initialized trials */
			LZe_backward(e, cur);
			return cur;
		}
		newlen = LZe_read_match_distances(e);
		if (newlen >= e->matchfinder->match_len_limit) {
			e->longest_match_found = newlen;
			LZe_backward(e, cur);
			return cur;
		}

		/* give final values to current trial */
		cur_trial = &e->trials[cur];
		prev_index = cur_trial->prev_index;
		cur_trial->state = e->trials[prev_index].state;

		for (i = 0; i < num_rep_distances; ++i)
			cur_trial->reps[i] = e->trials[prev_index].reps[i];

		if (prev_index == cur - 1) {
			if (cur_trial->dis == 0)
				cur_trial->state = St_set_short_rep(cur_trial->state);
			else
				cur_trial->state = St_set_char(cur_trial->state);
		} else {
			if (cur_trial->dis < num_rep_distances)
				cur_trial->state = St_set_rep(cur_trial->state);
			else
				cur_trial->state = St_set_match(cur_trial->state);
			mtf_reps(cur_trial->dis, cur_trial->reps);
		}

		pos_state = Mf_data_position(e->matchfinder) & pos_state_mask;
		prev_byte = Mf_peek(e->matchfinder, 1);
		cur_byte = Mf_peek(e->matchfinder, 0);
		match_byte = Mf_peek(e->matchfinder, cur_trial->reps[0] + 1);

		next_price = cur_trial->price +
			price0(e->bm_match[cur_trial->state][pos_state]);
		if (St_is_char(cur_trial->state))
			next_price += LZe_price_literal(e, prev_byte, cur_byte);
		else
			next_price += LZe_price_matched(e, prev_byte, cur_byte,
							match_byte);
		/* try last updates to next trial */
		next_trial = &e->trials[cur + 1];

		Tr_update(next_trial, next_price, -1, cur);	/* literal */

		match_price = cur_trial->price +
			price1(e->bm_match[cur_trial->state][pos_state]);
		rep_match_price = match_price + price1(e->bm_rep[cur_trial->state]);

		if (match_byte == cur_byte && next_trial->dis != 0)
			Tr_update(next_trial, rep_match_price +
				LZe_price_shortrep(e, cur_trial->state,
					pos_state), 0, cur);

		len_limit = MIN(MIN(max_num_trials - 1 - cur,
				Mf_available_bytes(e->matchfinder)),
				e->matchfinder->match_len_limit);
		if (len_limit < min_match_len) continue;

		for (rep = 0; rep < num_rep_distances; ++rep) {
			const int dis = cur_trial->reps[rep] + 1;
			int len = 0;
			const uint8_t *const data =
				Mf_ptr_to_current_pos(e->matchfinder);
			while (len < len_limit && data[len] == data[len - dis])
				++len;
			if (len >= min_match_len) {
				const int price = rep_match_price +
					LZe_price_rep(e, rep, cur_trial->state, pos_state);
				while (num_trials < cur + len)
					e->trials[++num_trials].price = infinite_price;
				for (; len >= min_match_len; --len)
					Tr_update(&e->trials[cur + len], price +
						Lee_price(&e->rep_len_encoder, len,
							pos_state), rep, cur);
			}
		}

		if (newlen <= len_limit &&
			(newlen > min_match_len ||
			(newlen == min_match_len &&
			e->match_distances[min_match_len] < modeled_distances))) {
			const int normal_match_price = match_price +
				price0(e->bm_rep[cur_trial->state]);
			int len;
			int dis = e->match_distances[min_match_len];
			int len_state = get_len_state(min_match_len);
			int dis_price = infinite_price;

			while (num_trials < cur + newlen)
				e->trials[++num_trials].price = infinite_price;

			if (dis < modeled_distances)
				Tr_update(&e->trials[cur + min_match_len],
					normal_match_price +
					e->dis_prices[len_state][dis] +
					Lee_price(&e->match_len_encoder,
						min_match_len, pos_state),
					dis + num_rep_distances, cur);

			for (len = min_match_len + 1; len <= newlen; ++len) {
				if (dis != e->match_distances[len] ||
					len_state < len_states - 1) {
					dis = e->match_distances[len];
					len_state = get_len_state(len);
					dis_price = LZe_price_dis(e, dis, len_state);
				}
				Tr_update(&e->trials[cur + len],
					normal_match_price + dis_price +
					Lee_price(&e->match_len_encoder, len, pos_state),
					dis + num_rep_distances, cur);
			}
		}
	}
}


static bool LZe_encode_member(struct LZ_encoder *const e)
{
	const int dis_price_count =
		(e->matchfinder->match_len_limit > 12) ? 512 : 2048;
	int dis_price_counter = 0;
	int ahead, i;
	int reps[num_rep_distances];
	State state = 0;
	for (i = 0; i < num_rep_distances; ++i) reps[i] = 0;

	if (!Mf_finished(e->matchfinder)) {	/* encode first byte */
		const uint8_t prev_byte = 0;
		const uint8_t cur_byte = Mf_peek(e->matchfinder, 0);
		Re_encode_bit(&e->renc, &e->bm_match[state][0], 0);
		LZe_encode_literal(e, prev_byte, cur_byte);
		CRC32_update_byte(&e->crc, cur_byte);
		Mf_longest_match_len(e->matchfinder, 0);
		Mf_move_pos(e->matchfinder);
	}

	while (!Mf_finished(e->matchfinder)) {
		if (dis_price_counter <= 0) {
			LZe_update_distance_prices(e);
			dis_price_counter = dis_price_count;
		}

		ahead = LZe_sequence_optimizer(e, reps, state);
		dis_price_counter -= ahead;

		for (i = 0; ahead > 0;) {
			const int pos_state =
				(Mf_data_position(e->matchfinder) - ahead) & pos_state_mask;
			const int dis = e->trials[i].dis;
			const int len = e->trials[i].price;

			bool bit = (dis < 0 && len == 1);
			Re_encode_bit(&e->renc, &e->bm_match[state][pos_state], !bit);
			if (bit) {	/* literal byte */
				const uint8_t prev_byte = Mf_peek(e->matchfinder, ahead + 1);
				const uint8_t cur_byte = Mf_peek(e->matchfinder, ahead);
				CRC32_update_byte(&e->crc, cur_byte);
				if (St_is_char(state))
					LZe_encode_literal(e, prev_byte, cur_byte);
				else {
					const uint8_t match_byte =
						Mf_peek(e->matchfinder, ahead + reps[0] + 1);
					LZe_encode_matched(e, prev_byte, cur_byte, match_byte);
				}
				state = St_set_char(state);
			} else {	/* match or repeated match */

				e->crc = crc32_block_endian0(e->crc,
						Mf_ptr_to_current_pos(e->matchfinder) - ahead,
						len, global_crc32_table);
				mtf_reps(dis, reps);
				bit = (dis < num_rep_distances);
				Re_encode_bit(&e->renc, &e->bm_rep[state], bit);
				if (bit) {		/* repeated match */
					bit = (dis == 0);
					Re_encode_bit(&e->renc, &e->bm_rep0[state], !bit);
					if (bit)
						Re_encode_bit(&e->renc, &e->bm_len[state][pos_state], len > 1);
					else {
						Re_encode_bit(&e->renc, &e->bm_rep1[state], dis > 1);
						if (dis > 1)
							Re_encode_bit(&e->renc, &e->bm_rep2[state], dis > 2);
					}
					if (len == 1)
						state = St_set_short_rep(state);
					else {
						Lee_encode(&e->rep_len_encoder,
							&e->renc, len, pos_state);
						state = St_set_rep(state);
					}
				} else {		/* match */
					LZe_encode_pair(e, dis - num_rep_distances, len, pos_state);
					state = St_set_match(state);
				}
			}
			ahead -= len;
			i += len;
		}
	}
	LZe_full_flush(e, state);
	return !e->renc.write_error;
}


struct Lzma_options {
	int dictionary_size;	/* 4KiB..512MiB */
	int match_len_limit;	/* 5..273 */
} encoder_options;


static int getnum(const char *const ptr, const int llimit, const int ulimit)
{
	long result;
	char *tail;
	errno = 0;
	result = strtol(ptr, &tail, 0);
	if (tail == ptr || errno)
		goto error;
	if (tail[0]) {
		int factor = (tail[1] == 'i') ? 1024 : 1000;
		int exponent = 0, i;

		switch (tail[0]) {
		case 'M':
			exponent = 2;
			break;
		case 'K':
			if (factor == 1024) {
				exponent = 1;
				break;
			}
			goto error;
		case 'k':
			if (factor == 1000) {
				exponent = 1;
				break;
			}
		default:
			goto error;
		}
		for (i = 0; i < exponent; ++i) {
			if (LONG_MAX / factor >= labs(result))
				result *= factor;
			else
				goto error;
		}
	}
	if (result >= llimit && result <= ulimit)
		return result;
  error:
	bb_error_msg_and_die("invalid number");
}


static int get_dict_size(const char *const arg)
{
	char *tail;
	long bits = strtol(arg, &tail, 0);
	if (bits >= min_dictionary_bits &&
		bits <= max_dictionary_bits && *tail == 0)
		return (1 << bits);
	return getnum(arg, min_dictionary_size, max_dictionary_size);
}


static IF_DESKTOP(long long) int FAST_FUNC pack_lzip(transformer_state_t *xstate UNUSED_PARAM)
{
	int retval = 0;
	File_header header;
	struct Matchfinder matchfinder;
	struct LZ_encoder * encoder;

	Fh_set_magic(header);
	if (!Fh_set_dictionary_size(header, encoder_options.dictionary_size) ||
		encoder_options.match_len_limit < min_match_len_limit ||
		encoder_options.match_len_limit > max_match_len)
		bb_error_msg_and_die("internal error");

	if (!Mf_init(&matchfinder, Fh_get_dictionary_size(header),
			encoder_options.match_len_limit)) {
		bb_error_msg(bb_msg_memory_exhausted);
		return -1;
	}
	Fh_set_dictionary_size(header, matchfinder.dictionary_size);

	encoder = malloc(sizeof(struct LZ_encoder));
	if (!encoder || !LZe_init(encoder, &matchfinder, header)) {
		bb_error_msg(bb_msg_memory_exhausted);
		retval = -1;
	} else {
		if (!LZe_encode_member(encoder)) {
			bb_perror_msg(bb_msg_write_error);
			retval = -1;
		}
		LZe_free(encoder);
	}
	free(encoder);
	Mf_free(&matchfinder);
	return retval;
}


int lzip_main(int argc, char **argv) MAIN_EXTERNALLY_VISIBLE;
int lzip_main(int argc UNUSED_PARAM, char **argv)
{
	/* Mapping from gzip/bzip2 style 1..9 compression modes
	   to the corresponding LZMA compression modes. */
	const struct Lzma_options option_mapping[] = {
		{1 << 20, 5},	/* -0 */
		{1 << 20, 5},	/* -1 */
		{3 << 19, 6},	/* -2 */
		{1 << 21, 8},	/* -3 */
		{3 << 20, 12},	/* -4 */
		{1 << 22, 20},	/* -5 */
		{1 << 23, 36},	/* -6 */
		{1 << 24, 68},	/* -7 */
		{3 << 23, 132},	/* -8 */
		{1 << 25, 273}	/* -9 */
	};
	int i;
	char *m_arg;
	char *s_arg;
	/* Must match bbunzip's constants OPT_STDOUT, OPT_FORCE! */
	uint32_t flags = getopt32(argv, "cfvqdt0123456789Fkm:s:", &m_arg, &s_arg);

	if (flags & 0x30) {	// -d and/or -t
#if ENABLE_LUNZIP		/* lunzip_main may not be visible... */
		return lunzip_main(argc, argv);
#else
		bb_error_msg("decompression is disabled");
		return 1;
#endif
	}
	flags >>= 6;

	encoder_options = option_mapping[6];	/* default = "-6" */

	for (i = 9; i >= 7; --i)
		if (flags & (1 << i))
			encoder_options = option_mapping[i];
	for (i = 0; i <= 6; ++i)
		if (flags & (1 << i))
			encoder_options = option_mapping[i];
	if (flags & (1 << 12))			/* -m */
		encoder_options.match_len_limit =
			getnum(m_arg, min_match_len_limit, max_match_len);
	if (flags & (1 << 13))			/* -s */
		encoder_options.dictionary_size = get_dict_size(s_arg);
	/* end process options */

	argv += optind;

	if (!global_crc32_table)
		global_crc32_table = crc32_filltable(NULL, 0);
	if (!dis_slots) {
		Dis_slots_init();
		Prob_prices_init();
	}

	return bbunpack(argv, pack_lzip, append_ext, "lz");
}
