#ifndef SPLIT_MAP_H
#define SPLIT_MAP_H

#include "frag_check.h"

typedef struct {
	int ref_start;
	int ref_end;
	int read_start;
	int read_end;
} hash_frag_aln;

typedef struct{
	int frag_num;
	int frag_max;
	hash_frag_aln *hfa;
} hash_frag;

typedef struct {
	int32_t offset;
	int32_t len;
}hash_blank_t;

typedef struct {
	int32_t n_blanks;
	hash_blank_t *blank;
}hash_line_t;

#define HASH_FRAG_START 1
#define HASH_FRAG_END 2
#define HASH_FRAG_SEED 3

//XXX
//64-slen or 32-slen
#define LAST_SHIFT(x, slen) ((x<<slen)>>slen)
extern const int8_t sc_mat[25];

int split_delete_map(uint32_t **res_cigar, int *res_len, uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, int64_t ref_offset, 
					int hash_len, uint32_t **hash_num, uint64_t ***hash_node, int key_len, int hash_size);
int split_insert_map(uint32_t **res_cigar, int *res_len, uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, int64_t ref_offset, 
					int hash_len, uint32_t **hash_num, uint64_t ***hash_node, int key_len, int hash_size);
#endif
