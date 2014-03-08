#ifndef FRAG_H
#define FRAG_H

#include <stdint.h>
#include "lsat_aln.h"
#include "bntseq.h"

/*
typedef struct {
	int frag_max;	//
	int frag_num;	//num of frag
	int last_len;	//0~2*PER_LEN-1 XXX
	int seed_max;
	int seed_num;	//num of seed
	int *chr;		
	int *srand;		//1:+, -1:-

	int *ref_begin;	//original exact pos, 1-base
	int *ref_end;
	int *read_begin;
	int *read_end;

	int *ex_ref_begin;	//extended exact pos, 1-base
	int *ex_ref_end;
	int *ex_read_begin;
	int *ex_read_end;

	int *per_n;		//PER_LEN num
	
	int *flag;		//COVERED or UNCOVERED
	int *b_f;		//boundary frag

	//seed msg
	int **seed_i;	//index of aln_msg
	int **seed_aln_i;//index of *aln_msg
} frag_msg;
*/

typedef struct {
	int chr;
	int srand;

	int ref_begin;	//1-based   //XXX are they useful?
	int ref_end;
	int read_begin;	//1-based
	int read_end;
	
	int ex_ref_begin;	//1-based
	int ex_ref_end;
	int ex_read_begin;	//1-based
	int ex_read_end;
	
	int64_t cigar_ref_start;	//frag's cigar start of ref, 1-based
	int64_t cigar_ref_end;		
	int64_t cigar_read_start;	//frag's cigar start of read, 1-based
	int64_t cigar_read_end;
	uint32_t *cigar;		//frag's cigar
	int cigar_len;			//frag's cigar length
	int cigar_max;			//size of cigar[]
	int edit_dis;			//frag's edit-dis
	int len_dif;			//length difference between ref and read. eg, ref=101, read=100, then len_dif = 101-100 = 1.

	int per_n;
	int flag;
	int b_f;
	
	int seed_num;	//seed number of per frag
	int *seed_i;
	int *seed_aln_i;
} frag_aln_msg;

typedef struct {
	int frag_max;
	int frag_num;
	int last_len;
	int seed_num;	//whole number
	int per_seed_max;
	frag_aln_msg *fa_msg;
} frag_msg;


typedef struct {
	uint64_t offset;	//1-based
	int chr;
	int nsrand;			//1:'+' -1:'-'
	uint32_t *cigar;
	int c_m;
	int cigar_len;
} aln_res;

extern const int8_t sc_mat[25];
extern const int8_t bwasw_sc_mat[25];
extern char nst_nt4_table[256];
frag_msg *frag_init_msg(int frag_max);
void frag_free_msg(frag_msg *f_msg);
int frag_set_msg(aln_msg *a_msg, int seed_i, int aln_i, int FLAG, frag_msg *f_msg, int frag_i, int seed_len);//FLAG 0: start/1:end / 2:seed
int frag_check(char *read_name, bntseq_t *bns, uint8_t *pac, const char *read_prefix, char *read_seq, int read_len, int seed_all, frag_msg *f_msg, aln_msg *a_msg, uint32_t **hash_num, uint64_t ***hash_node, int seed_len, int last_len);
void printcigar(uint32_t *cigar, int cigar_len);

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

#endif
