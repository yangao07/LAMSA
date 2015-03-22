#ifndef BWT_ALN_H
#define BWT_ALN_H
#include "lsat_aln.h"
#include "bwt.h"
#include "bntseq.h"


typedef struct {
    int ref_id;
    uint8_t is_rev;
    uint64_t ref_pos;
    line_node from;
    int score;
} loc_t;

typedef struct {
    int pos;
    loc_t *loc;
    int n, m;
} bwt_seed_t;

typedef struct {
    int read_pos;
    uint64_t ref_pos;
} bwt_bound;

void bwt_aln_remain(aln_reg *a_reg, aln_res *a_res, bwt_t *bwt, bntseq_t *bns, uint8_t *pac, char *read_seq, lsat_aln_per_para *APP, lsat_aln_para *AP);

#endif
