#ifndef FRAG_H
#define FRAG_H

#include <stdint.h>
#include "lsat_aln.h"
#include "bntseq.h"

typedef struct {
	int chr;
	int srand;

	int64_t cigar_ref_start;	//frag's cigar start of ref, 1-based    XXX
	int64_t cigar_ref_end;		                                    //  XXX
	int64_t cigar_read_start;	//frag's cigar start of read, 1-based   XXX
	int64_t cigar_read_end;                                         //  XXX 
	int32_t *cigar;		//frag's cigar
	int cigar_len;			//frag's cigar length
	int cigar_max;			//size of cigar[]
	//int edit_dis;			//frag's edit-dis
	int len_dif;			//length difference between ref and read. eg, ref=101, read=100, then len_dif = 101-100 = 1.

	int per_n;
	int flag;
	int b_f;
	
	int seed_num;	//seed number of per frag
	int *seed_i;
	int *seed_aln_i;

    int next_trigger_n;
    int pre_trigger_n;
    int trigger_m;
    int *trigger;
} frag_aln_msg;

typedef struct {
	int frag_max;
	int frag_num;
	//int last_len; (APP.last_len)
	//int seed_num;	//whole number (APP.seed_out)
	//int seed_all;	//all the seeds from read (APP.seed_all)
	int per_seed_max;
	frag_aln_msg *fa_msg;
	int frag_left_bound;	//read_id of previous line's last seed, if NO pre line, this will be 0
	int frag_right_bound;	//read_id of next line's first seed, if NO next line, this will be n_seed+1
} frag_msg;

typedef struct {
	uint64_t offset;	//1-based
	int chr;
	int nsrand;			//1:'+' 0:'-'
	int32_t *cigar;
	int c_m;
	int cigar_len;
} line_aln_res;

typedef struct {
	int res_m;
	int cur_res_n;
	line_aln_res *la;	
    int read_len;
} aln_res;

extern const int8_t sc_mat[25];
extern const int8_t bwasw_sc_mat[25];
extern char nst_nt4_table[256];
void _push_cigar(int32_t **cigar, int *cigar_len, int *cigar_m, int32_t *_cigar, int _cigar_len);
void _push_cigar1(int32_t **cigar, int *cigar_len, int *cigar_m, int32_t _cigar);
void frag_init_msg(frag_msg *f_msg, int frag_max);
void frag_free_msg(frag_msg *f_msg, int line_num);
int frag_set_msg(aln_msg *a_msg, int seed_i, int aln_i, int FLAG, frag_msg *f_msg, int frag_i);//FLAG 0: start/1:end / 2:seed
int frag_trigger_set(frag_dp_node f_node, frag_msg *f_msg, int frag_i);
int frag_copy_msg(frag_msg *ff_msg, frag_msg *tf_msg);
int frag_check(aln_msg *a_msg, frag_msg **f_msg,
               bntseq_t *bns, uint8_t *pac, const char *read_prefix,
               char *read_seq,
               lsat_aln_per_para *APP, lsat_aln_para *AP,
               int line_n, int *line_tri,
               uint32_t **hash_num, uint64_t ***hash_node);

void printcigar(FILE *outp, int32_t *cigar, int cigar_len);

#define MAXOFTWO(a, b) ((a) > (b) ? (a) : (b))
#define MINOFTWO(a, b) ((a) < (b) ? (a) : (b))
#define MAX_FRAG 2048 
#define MAX_LINE_LEN 65536

#define GOOD_ALIGN 0
#define BAD_ALIGN 1

#define RIGHT 1
#define LEFT  0

#define COVERED 1
#define UNCOVERED 0

#define MAX_E 3

#define _MAX_K 30

#endif
