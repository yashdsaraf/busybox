/*
 * lunzip implementation for busybox
 *
 * Copyright (C) 2012-2015 Antonio Diaz Diaz.
 *
 * Licensed under GPLv2 or later, see file LICENSE in this source tree.
 */

#include "libbb.h"
#include "bb_archive.h"
#include "lzip.h"


enum { rd_buffer_size = 16384 };

struct Range_decoder {
	unsigned long long partial_member_pos;
	uint8_t *buffer;	/* input buffer */
	int pos;		/* current pos in buffer */
	int stream_pos;		/* when reached, a new block must be read */
	uint32_t code;
	uint32_t range;
	int infd;		/* input file descriptor */
	bool at_stream_end;
};


static bool Rd_read_block(struct Range_decoder *const rdec)
{
	if (!rdec->at_stream_end) {
		rdec->stream_pos =
			full_read(rdec->infd, rdec->buffer, rd_buffer_size);
		rdec->at_stream_end = (rdec->stream_pos < rd_buffer_size);
		rdec->partial_member_pos += rdec->pos;
		rdec->pos = 0;
	}
	return rdec->pos < rdec->stream_pos;
}


static bool Rd_init(struct Range_decoder *const rdec, const int ifd,
			const bool preload_magic)
{
	rdec->partial_member_pos = (preload_magic ? 4 : 0);
	rdec->buffer = (uint8_t *) malloc(rd_buffer_size);
	if (!rdec->buffer) return false;
	rdec->pos = 0;
	rdec->stream_pos = 0;
	rdec->code = 0;
	rdec->range = 0xFFFFFFFFU;
	rdec->infd = ifd;
	rdec->at_stream_end = false;
	return true;
}

static inline void Rd_free(struct Range_decoder *const rdec)
{
	free(rdec->buffer);
}


static inline bool Rd_finished(struct Range_decoder *const rdec)
{
	return rdec->pos >= rdec->stream_pos && !Rd_read_block(rdec);
}

static inline unsigned long long
Rd_member_position(const struct Range_decoder *const rdec)
{
	return rdec->partial_member_pos + rdec->pos;
}

static inline void Rd_reset_member_position(struct Range_decoder *const rdec)
{
	rdec->partial_member_pos = -rdec->pos;
}

static uint8_t Rd_get_byte(struct Range_decoder *const rdec)
{
	if (Rd_finished(rdec)) return 0;
	return rdec->buffer[rdec->pos++];
}

static void Rd_load(struct Range_decoder *const rdec)
{
	int i;
	rdec->code = 0;
	for (i = 0; i < 5; ++i)
		rdec->code = (rdec->code << 8) | Rd_get_byte(rdec);
	rdec->range = 0xFFFFFFFFU;
}

static void Rd_normalize(struct Range_decoder *const rdec)
{
	if (rdec->range <= 0x00FFFFFFU) {
		rdec->range <<= 8;
		rdec->code = (rdec->code << 8) | Rd_get_byte(rdec);
	}
}

static int Rd_decode(struct Range_decoder *const rdec, const int num_bits)
{
	int symbol = 0;
	int i;
	for (i = num_bits; i > 0; --i) {
		uint32_t mask;
		Rd_normalize(rdec);
		rdec->range >>= 1;
		/* symbol <<= 1; */
		/* if(rdec->code >= rdec->range) { rdec->code -= rdec->range; symbol |= 1; } */
		mask = 0U - (rdec->code < rdec->range);
		rdec->code -= rdec->range;
		rdec->code += rdec->range & mask;
		symbol = (symbol << 1) + (mask + 1);
	}
	return symbol;
}

static int Rd_decode_bit(struct Range_decoder *const rdec,
			Bit_model * const probability)
{
	uint32_t bound;
	Rd_normalize(rdec);
	bound = (rdec->range >> bit_model_total_bits) * *probability;
	if (rdec->code < bound) {
		rdec->range = bound;
		*probability += (bit_model_total - *probability) >> bit_model_move_bits;
		return 0;
	} else {
		rdec->range -= bound;
		rdec->code -= bound;
		*probability -= *probability >> bit_model_move_bits;
		return 1;
	}
}

static int Rd_decode_tree(struct Range_decoder *const rdec,
				Bit_model bm[], const int num_bits)
{
	int symbol = 1;
	int i;
	for (i = num_bits; i > 0; --i)
		symbol = (symbol << 1) | Rd_decode_bit(rdec, &bm[symbol]);
	return symbol - (1 << num_bits);
}

static int Rd_decode_tree_reversed(struct Range_decoder *const rdec,
					Bit_model bm[], const int num_bits)
{
	int model = 1;
	int symbol = 0;
	int i;
	for (i = 0; i < num_bits; ++i) {
		const bool bit = Rd_decode_bit(rdec, &bm[model]);
		model <<= 1;
		if (bit) { ++model; symbol |= (1 << i); }
	}
	return symbol;
}

static int Rd_decode_matched(struct Range_decoder *const rdec,
				Bit_model bm[], int match_byte)
{
	Bit_model *const bm1 = bm + 0x100;
	int symbol = 1;
	while(symbol < 0x100) {
		int match_bit, bit;
		match_byte <<= 1;
		match_bit = match_byte & 0x100;
		bit = Rd_decode_bit(rdec, &bm1[match_bit+symbol]);
		symbol = (symbol << 1) | bit;
		if (match_bit != bit << 8) {
			while(symbol < 0x100)
				symbol = (symbol << 1) | Rd_decode_bit(rdec, &bm[symbol]);
			break;
		}
	}
	return symbol & 0xFF;
}

static int Rd_decode_len(struct Range_decoder *const rdec,
			struct Len_model * const lm,
			const int pos_state)
{
	if (Rd_decode_bit(rdec, &lm->choice1) == 0)
		return Rd_decode_tree(rdec, lm->bm_low[pos_state], len_low_bits);
	if (Rd_decode_bit(rdec, &lm->choice2) == 0)
		return len_low_symbols +
			Rd_decode_tree(rdec, lm->bm_mid[pos_state], len_mid_bits);
	return len_low_symbols + len_mid_symbols +
		Rd_decode_tree(rdec, lm->bm_high, len_high_bits);
}


struct LZ_decoder {
	unsigned long long partial_data_pos;
	struct Range_decoder *rdec;
	int buffer_size;	/* buffer_size == dictionary_size */
	uint8_t *buffer;	/* output buffer */
	int pos;		/* current pos in buffer */
	int stream_pos;		/* first byte not yet written to file */
	uint32_t crc;
	int outfd;		/* output file descriptor */

	Bit_model bm_literal[1 << literal_context_bits][0x300];
	Bit_model bm_match[states][pos_states];
	Bit_model bm_rep[states];
	Bit_model bm_rep0[states];
	Bit_model bm_rep1[states];
	Bit_model bm_rep2[states];
	Bit_model bm_len[states][pos_states];
	Bit_model bm_dis_slot[len_states][1 << dis_slot_bits];
	Bit_model bm_dis[modeled_distances-end_dis_model];
	Bit_model bm_align[dis_align_size];

	struct Len_model match_len_model;
	struct Len_model rep_len_model;
	bool write_error;
};

static void LZd_flush_data(struct LZ_decoder *const d)
{
	if (d->pos > d->stream_pos) {
		const int size = d->pos - d->stream_pos;
		d->crc = crc32_block_endian0(d->crc, d->buffer + d->stream_pos,
						size, global_crc32_table);
		if (d->outfd >= 0 && full_write(d->outfd,
				d->buffer + d->stream_pos, size) != size)
			d->write_error = true;
		if (d->pos >= d->buffer_size) {
			d->partial_data_pos += d->pos;
			d->pos = 0;
		}
		d->stream_pos = d->pos;
	}
}

static uint8_t LZd_get_prev_byte(const struct LZ_decoder *const d)
{
	const int i = ((d->pos > 0) ? d->pos : d->buffer_size) - 1;
	return d->buffer[i];
}

static uint8_t LZd_get_byte(const struct LZ_decoder *const d,
				const int distance)
{
	int i = d->pos - distance - 1;
	if (i < 0) i += d->buffer_size;
	return d->buffer[i];
}

static void LZd_put_byte(struct LZ_decoder *const d, const uint8_t b)
{
	d->buffer[d->pos] = b;
	if (++d->pos >= d->buffer_size) LZd_flush_data(d);
}

static void LZd_copy_block(struct LZ_decoder *const d,
				const int distance, int len)
{
	int i = d->pos - distance - 1;
	if (i < 0) i += d->buffer_size;
	if (len < d->buffer_size - MAX(d->pos, i) && len <= abs(d->pos - i)) {
		/* no wrap, no overlap */
		memcpy(d->buffer + d->pos, d->buffer + i, len);
		d->pos += len;
	} else for (; len > 0; --len) {
		d->buffer[d->pos] = d->buffer[i];
		if (++d->pos >= d->buffer_size) LZd_flush_data(d);
		if (++i >= d->buffer_size) i = 0;
	}
}

static bool LZd_init(struct LZ_decoder *const d,
			struct Range_decoder *const rde,
			const int dict_size, const int ofd)
{
	d->partial_data_pos = 0;
	d->rdec = rde;
	d->buffer_size = dict_size;
	d->buffer = (uint8_t *) malloc(d->buffer_size);
	if (!d->buffer) return false;
	d->pos = 0;
	d->stream_pos = 0;
	d->crc = 0xFFFFFFFFU;
	d->outfd = ofd;
	Bm_array_init(&d->bm_literal[0][0], lz_num_models);
	Lm_init(&d->match_len_model);
	Lm_init(&d->rep_len_model);
	d->write_error = false;
	d->buffer[d->buffer_size - 1] = 0;	/* prev_byte of first byte */
	return true;
}

static inline void LZd_free(struct LZ_decoder *const d)
{
	free(d->buffer);
}

static inline unsigned LZd_crc(const struct LZ_decoder *const d)
{
	return d->crc ^ 0xFFFFFFFFU;
}

static inline unsigned long long
LZd_data_position(const struct LZ_decoder *const d)
{
	return d->partial_data_pos + d->pos;
}


static bool LZd_verify_trailer(struct LZ_decoder *const d)
{
	File_trailer trailer;
	int i = 0;
	while (i < Ft_size)
		trailer[i++] = Rd_get_byte(d->rdec);

	return (Ft_get_data_crc(trailer) == LZd_crc(d) &&
		Ft_get_data_size(trailer) == LZd_data_position(d) &&
		Ft_get_member_size(trailer) == Rd_member_position(d->rdec));
}


/* Return value: -1 = write error, 0 = OK, 1 = data error. */
static int LZd_decode_member(struct LZ_decoder *const d)
{
	struct Range_decoder * const rdec = d->rdec;
	unsigned rep0 = 0;	/* rep[0-3] latest four distances */
	unsigned rep1 = 0;	/* used for efficient coding of */
	unsigned rep2 = 0;	/* repeated distances */
	unsigned rep3 = 0;
	State state = 0;

	Rd_load(rdec);
	while (!Rd_finished(rdec)) {
		const int pos_state = LZd_data_position(d) & pos_state_mask;
		if (Rd_decode_bit(rdec, &d->bm_match[state][pos_state]) == 0) {
			const uint8_t prev_byte = LZd_get_prev_byte(d);
			if (St_is_char(state)) {
				state -= (state < 4) ? state : 3;
				LZd_put_byte(d, Rd_decode_tree(rdec,
						d->bm_literal[get_lit_state(prev_byte)], 8));
			} else {
				state -= (state < 10) ? 3 : 6;
				LZd_put_byte(d, Rd_decode_matched(rdec,
						d->bm_literal[get_lit_state(prev_byte)],
						LZd_get_byte(d, rep0)));
			}
		} else {
			int len;
			if (Rd_decode_bit(rdec, &d->bm_rep[state]) != 0) {
				if (Rd_decode_bit(rdec, &d->bm_rep0[state]) != 0) {
					unsigned distance;
					if (Rd_decode_bit(rdec, &d->bm_rep1[state]) == 0)
						distance = rep1;
					else {
						if (Rd_decode_bit(rdec, &d->bm_rep2[state]) == 0)
							distance = rep2;
						else {
							distance = rep3;
							rep3 = rep2;
						}
						rep2 = rep1;
					}
					rep1 = rep0;
					rep0 = distance;
				} else {
					if (Rd_decode_bit(rdec, &d->bm_len[state][pos_state]) == 0) {
						state = St_set_short_rep(state);
						LZd_put_byte(d, LZd_get_byte(d, rep0));
						continue;
					}
				}
				state = St_set_rep(state);
				len = min_match_len + Rd_decode_len(rdec, &d->rep_len_model, pos_state);
			} else {
				int dis_slot;
				const unsigned rep0_saved = rep0;
				len = min_match_len + Rd_decode_len(rdec, &d->match_len_model, pos_state);
				dis_slot = Rd_decode_tree(rdec, d->bm_dis_slot[get_len_state(len)], 6);
				if (dis_slot < start_dis_model) rep0 = dis_slot;
				else {
					const int direct_bits = (dis_slot >> 1) - 1;
					rep0 = (2 | (dis_slot & 1)) << direct_bits;
					if (dis_slot < end_dis_model)
						rep0 += Rd_decode_tree_reversed(rdec,
							d->bm_dis + rep0 - dis_slot - 1, direct_bits);
					else {
						rep0 +=	Rd_decode(rdec, direct_bits - dis_align_bits) << dis_align_bits;
						rep0 += Rd_decode_tree_reversed(rdec, d->bm_align, dis_align_bits);
						if (rep0 == 0xFFFFFFFFU) {	/* Marker found */
							rep0 = rep0_saved;
							Rd_normalize(rdec);
							LZd_flush_data(d);
							if (d->write_error) return -1;
							if (len == min_match_len &&	/* End Of Stream marker */
								LZd_verify_trailer(d))
								return 0;
							if (len == min_match_len + 1) {	/* Sync Flush marker */
								Rd_load(rdec);
								continue;
							}
							return 1;
						}
					}
				}
				rep3 = rep2;
				rep2 = rep1;
				rep1 = rep0_saved;
				state = St_set_match(state);
				if (rep0 >= (unsigned) d->buffer_size ||
					rep0 >= LZd_data_position(d)) {
					LZd_flush_data(d);
					return 1;
				}
			}
			LZd_copy_block(d, rep0, len);
		}
	}
	LZd_flush_data(d);
	return 1;
}


IF_DESKTOP(long long) int FAST_FUNC
unpack_lz_stream(transformer_state_t *xstate)
{
	IF_DESKTOP(long long) int total = 0;
	struct Range_decoder rdec;
	bool first_member;
	const bool preload_magic = (xstate->check_signature == 0);

	if (!global_crc32_table)
		global_crc32_table = crc32_filltable(NULL, 0);

	if (!Rd_init(&rdec, xstate->src_fd, preload_magic))
		return -1;

	for (first_member = true;; first_member = false) {
		int tmp = 0;
		File_header header;
		struct LZ_decoder * decoder;

		if (first_member && preload_magic) {
			Fh_set_magic(header);
			tmp = 4;
		} else {
			Rd_reset_member_position(&rdec);
		}
		while (tmp < Fh_size)
			header[tmp++] = Rd_get_byte(&rdec);
		if (Rd_finished(&rdec)) {	/* End Of File */
			if (first_member) {
				bb_error_msg(bb_msg_read_error);
				total = -1;
			}
			break;
		}
		tmp = Fh_get_dictionary_size(header);
		if (!Fh_verify_magic(header) || tmp < min_dictionary_size ||
			tmp > max_dictionary_size) {
			if (!first_member)
				break;		/* trailing garbage */
			bb_error_msg("invalid magic");
			total = -1;
			break;
		}

		decoder = malloc(sizeof(struct LZ_decoder));
		if (!decoder || !LZd_init(decoder, &rdec, tmp, xstate->dst_fd)) {
			bb_error_msg(bb_msg_memory_exhausted);
			free(decoder);
			total = -1;
			break;
		}
		tmp = LZd_decode_member(decoder);
		IF_DESKTOP(total += Rd_member_position(&rdec);)
		LZd_free(decoder);
		free(decoder);
		if (tmp != 0) {
			if (tmp < 0)
				bb_perror_msg(bb_msg_write_error);
			else
				bb_error_msg("corrupted data");
			total = -1;
			break;
		}
	}
	Rd_free(&rdec);
	return total;
}
