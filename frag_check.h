#ifndef FRAG_H
#define FRAG_H

#include <stdint.h>

typedef struct {
	int frag_max;
	int frag_num;
	int last_len;	//0-2*PER_LEN-1
	int read_num;
	int *chr;
	int *srand;		//1:+, -1:-

	int *read_1;	//soap2-dp	won't be used!
	int *read_2;
	int *ref_1;
	int *ref_2;

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
} frag_msg;

frag_msg *frag_init_msg(int frag_max);
void frag_free_msg(frag_msg *f_msg);

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
