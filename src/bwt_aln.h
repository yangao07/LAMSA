#ifndef BWT_ALN_H
#define BWT_ALN_H
#include "lamsa_aln.h"
#include "frag_check.h"
#include "bwt.h"
#include "bntseq.h"


typedef struct {
    int ref_id;
    uint8_t is_rev;
    ref_pos_t ref_pos;
    line_node from;
    int score, NM, node_n;
    int track_flag;
} loc_t;

typedef struct {
    int pos;
    loc_t *loc;
    int n, m;
} bwt_seed_t;

typedef struct {
    int read_pos;
    ref_pos_t ref_pos;
} bwt_bound;

void bwt_aln_remain(aln_reg *a_reg, aln_res *re_res, bwt_t *bwt, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, uint8_t **read_rbseq, lamsa_aln_para *AP, kseq_t *seqs);

#endif
