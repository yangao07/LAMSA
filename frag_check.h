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
	int ref_begin;
	int ref_end;
	int read_begin;
	int read_end;
	int ex_ref_begin;
	int ex_ref_end;
	int ex_read_begin;
	int ex_read_end;
	int per_n;
	int flag;
	int b_f;
	int seed_num;	//seed number of per frag
	int *seed_i;
	int *seed_aln_i;
} frag_aln_msg;

typedef struct {
	int frag_max;	//macro	XXX
	int frag_num;
	int last_len;
	int seed_num;	//whole number
	int per_seed_max;
	frag_aln_msg *fa_msg;
} frag_msg;


extern char nst_nt4_table[256];
frag_msg *frag_init_msg(int frag_max);
void frag_free_msg(frag_msg *f_msg);
int frag_set_msg(aln_msg *a_msg, int seed_i, int aln_i, int FLAG, frag_msg *f_msg, int frag_i, int seed_len);//FLAG 0: start/1:end / 2:seed
int frag_check(bntseq_t *bns, int8_t *pac, char *read_prefix, char *read_seq, frag_msg *f_msg, int seed_len);

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
