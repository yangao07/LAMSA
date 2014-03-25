#ifndef SPLIT_MAP_H
#define SPLIT_MAP_H
//#define __DEBUG__
//#define __NEW__

#include "lsat_aln.h"
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

typedef struct {
	int32_t n_line;
	int32_t *line_i;
	int32_t ref_cover_len;	//whole covered length
	int32_t read_cover_len;
}hash_line;

typedef struct {
	int32_t node_i;

	hash_line line;
	int32_t ref_blank;
	int32_t read_blank;

	int32_t mis_match;	//XXX

	int32_t from_i;	//previous node_i
} hash_DP_node;

typedef struct {
	line_node from;
	int32_t ref_i;
	//	int32_t read_i;
	int32_t offset;		//hash_pos[]
	int32_t score;
	int32_t node_n;
	uint8_t match_flag;
	int dp_flag;
} hash_dp_node;

#define MIN_FLAG 1
#define MULTI_FLAG 2

#define HASH_FRAG_START 1
#define HASH_FRAG_END 2
#define HASH_FRAG_SEED 3

#define HASH_SV_PEN 10

//XXX
//64-slen or 32-slen
#define LAST_SHIFT(x, slen) ((x<<slen)>>slen)
extern const int8_t sc_mat[25];
#define HASH_LEN 10
#define HASH_STEP 1

int split_delete_map(uint32_t **res_cigar, int *res_len, uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, int64_t ref_offset, 
					int hash_len, int hash_step, uint32_t **hash_num, uint64_t ***hash_node, int key_len, int hash_size);
int split_insert_map(uint32_t **res_cigar, int *res_len, uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, int64_t ref_offset, 
					int hash_len, int hash_step, uint32_t **hash_num, uint64_t ***hash_node, int key_len, int hash_size);
#endif
