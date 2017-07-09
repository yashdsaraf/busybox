/*  Lzip - LZMA lossless data compressor
    Copyright (C) 2008-2016 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

typedef int State;

enum { states = 12 };

static inline bool St_is_char(const State st) { return st < 7; }

static inline State St_set_char(const State st)
{
	static const State next[states] = { 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 4, 5 };
	return next[st];
}

static inline State St_set_match(const State st)
{
	return ((st < 7) ? 7 : 10);
}

static inline State St_set_rep(const State st)
{
	return ((st < 7) ? 8 : 11);
}

static inline State St_set_short_rep(const State st)
{
	return ((st < 7) ? 9 : 11);
}


enum {
	min_dictionary_bits = 12,
	min_dictionary_size = 1 << min_dictionary_bits,	/* >= modeled_distances */
	max_dictionary_bits = 29,
	max_dictionary_size = 1 << max_dictionary_bits,
	literal_context_bits = 3,
	pos_state_bits = 2,
	pos_states = 1 << pos_state_bits,
	pos_state_mask = pos_states - 1,

	len_states = 4,
	dis_slot_bits = 6,
	start_dis_model = 4,
	end_dis_model = 14,
	modeled_distances = 1 << (end_dis_model / 2),	/* 128 */
	dis_align_bits = 4,
	dis_align_size = 1 << dis_align_bits,

	len_low_bits = 3,
	len_mid_bits = 3,
	len_high_bits = 8,
	len_low_symbols = 1 << len_low_bits,
	len_mid_symbols = 1 << len_mid_bits,
	len_high_symbols = 1 << len_high_bits,
	max_len_symbols = len_low_symbols + len_mid_symbols + len_high_symbols,

	min_match_len = 2,					/* must be 2 */
	max_match_len = min_match_len + max_len_symbols - 1,	/* 273 */
	min_match_len_limit = 5,

	lz_num_models =
		((1 << literal_context_bits) * 0x300) +
		(2 * states * pos_states) +
		(4 * states) +
		(len_states * (1 << dis_slot_bits)) +
		(modeled_distances - end_dis_model) +
		dis_align_size,
};

static inline int get_len_state(const int len)
{
	return MIN(len - min_match_len, len_states - 1);
}

static inline int get_lit_state(const uint8_t prev_byte)
{
	return (prev_byte >> (8 - literal_context_bits));
}


enum { bit_model_move_bits = 5,
	bit_model_total_bits = 11,
	bit_model_total = 1 << bit_model_total_bits
};

typedef int Bit_model;

static inline void Bm_init(Bit_model * const probability)
{
	*probability = bit_model_total / 2;
}

static inline void Bm_array_init(Bit_model * const p, const int size)
{
	int i = 0;
	while (i < size)
		p[i++] = bit_model_total / 2;
}

struct Len_model {
	Bit_model choice1;
	Bit_model choice2;
	Bit_model bm_low[pos_states][len_low_symbols];
	Bit_model bm_mid[pos_states][len_mid_symbols];
	Bit_model bm_high[len_high_symbols];
};

static inline void Lm_init(struct Len_model * const lm)
{
	Bm_init(&lm->choice1);
	Bm_init(&lm->choice2);
	Bm_array_init(lm->bm_low[0], pos_states * len_low_symbols);
	Bm_array_init(lm->bm_mid[0], pos_states * len_mid_symbols);
	Bm_array_init(lm->bm_high, len_high_symbols);
}


static inline int real_bits(unsigned value)
{
	int bits = 0;
	while(value > 0) { value >>= 1; ++bits; }
	return bits;
}


static const uint8_t magic_string[4] = { 0x4C, 0x5A, 0x49, 0x50 }; /* "LZIP" */

typedef uint8_t File_header[6];		/* 0-3 magic bytes */
					/*   4 version */
					/*   5 coded_dict_size */
enum { Fh_size = 6 };

static inline void Fh_set_magic(File_header data)
{
	memcpy(data, magic_string, 4);
	data[4] = 1;
}

static inline bool Fh_verify_magic(const File_header data)
{
	return (memcmp(data, magic_string, 4) == 0 && data[4] == 1);
}

static inline unsigned Fh_get_dictionary_size(const File_header data)
{
	unsigned sz = (1 << (data[5] & 0x1F));
	if (sz > min_dictionary_size)
		sz -= (sz / 16) * ((data[5] >> 5) & 7);
	return sz;
}

static inline bool Fh_set_dictionary_size(File_header data, const unsigned sz)
{
	if (sz < min_dictionary_size || sz > max_dictionary_size) return false;
	data[5] = real_bits(sz - 1);
	if (sz > min_dictionary_size) {
		const unsigned base_size = 1 << data[5];
		const unsigned fraction = base_size / 16;
		unsigned i;
		for (i = 7; i >= 1; --i)
			if (base_size - (i * fraction) >= sz) {
				data[5] |= (i << 5);
				break;
			}
	}
	return true;
}


typedef uint8_t File_trailer[20];
			/*  0-3  CRC32 of the uncompressed data */
			/*  4-11 size of the uncompressed data */
			/* 12-19 member size including header and trailer */

enum { Ft_size = 20 };

static inline unsigned Ft_get_data_crc(const File_trailer data)
{
	unsigned tmp = 0;
	int i;
	for (i = 3; i >= 0; --i) {
		tmp <<= 8;
		tmp += data[i];
	}
	return tmp;
}

static inline void Ft_set_data_crc(File_trailer data, unsigned crc)
{
	int i;
	for (i = 0; i <= 3; ++i) {
		data[i] = (uint8_t)crc;
		crc >>= 8;
	}
}

static inline unsigned long long Ft_get_data_size(const File_trailer data)
{
	unsigned long long tmp = 0;
	int i;
	for (i = 11; i >= 4; --i) {
		tmp <<= 8;
		tmp += data[i];
	}
	return tmp;
}

static inline void Ft_set_data_size(File_trailer data, unsigned long long sz)
{
	int i;
	for (i = 4; i <= 11; ++i) {
		data[i] = (uint8_t)sz;
		sz >>= 8;
	}
}

static inline unsigned long long Ft_get_member_size(const File_trailer data)
{
	unsigned long long tmp = 0;
	int i;
	for (i = 19; i >= 12; --i) {
		tmp <<= 8;
		tmp += data[i];
	}
	return tmp;
}

static inline void Ft_set_member_size(File_trailer data, unsigned long long sz)
{
	int i;
	for (i = 12; i <= 19; ++i) {
		data[i] = (uint8_t)sz;
		sz >>= 8;
	}
}
