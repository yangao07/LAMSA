#ifndef FRAG_H
#define FRAG_H


#include <stdint.h>
#include "lamsa_aln.h"
#include "gem_parse.h"
#include "bntseq.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
	int chr;
	int strand;

	uint64_t cigar_ref_start;	//frag's cigar start of ref, 1-based
	uint64_t cigar_ref_end;		                                    
	int cigar_read_start;	//frag's cigar start of read, 1-based   
	int cigar_read_end;                                         
	cigar32_t *cigar;		//frag's cigar
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
} frag_aln_msg;

typedef struct {
	int frag_max;
	int frag_num;

	int per_seed_max;
	frag_aln_msg *fa_msg;
    int line_score;
	int frag_left_bound;  //read_id of previous line's last seed, if NO pre line, this will be 0
	int frag_right_bound; //read_id of next line's first seed, if NO next line, this will be n_seed+1
} frag_msg;

typedef struct {
    uint64_t offset;	//1-based
    int chr;
    int nstrand;			//1:'+' 0:'-'
    cigar32_t *cigar;
    int c_m; int cigar_len;
    uint64_t refend; int readend; // for merge_cigar

    int score;          // alignment score
    int NM;             // edit distance
    //cigar32_t *MD;    // mis-match postion string XXX

    int reg_beg, reg_end; // set in push_reg_res()
} res_t;

typedef struct {
    int line_score;
    line_node merg_msg;    // after filtering:
                           // 1,x: keep
                           // 0,x: dump
    res_t **XA;
    int XA_n, XA_m;        // per_max_multi

    int res_m, cur_res_n;
    res_t *res;
    int tol_score, tol_NM; // tol score of all res, including split penalty 
    uint8_t mapQ;
} line_aln_res;

typedef struct {
	int l_m, l_n;
	line_aln_res *la;	

    int read_len;
    float cov_f;
} aln_res;
// cigar operatios
//pop_mode: 0 pop_n -> cigar_n
//          1 pop_n -> base n
int readInCigar(cigar32_t *cigar, int cigar_len);
int refInCigar(cigar32_t *cigar, int cigar_len);
static inline int _pop_cigar(cigar32_t **cigar, int *cigar_n, int pop_mode, int pop_n) {
    if (pop_mode == 0) {
        (*cigar_n) -= pop_n;
        if (*cigar_n < 0) *cigar_n = 0;
    } else if (pop_mode == 1) {
        int i;
        int pop_l = pop_n;

        i = *cigar_n - 1;
        while (pop_l != 0)
        {
            if ((((*cigar)[i]) & 0xf) == CMATCH || (((*cigar)[i]) & 0xf) == CINS) {
                if ((((*cigar)[i]) >> 4) < pop_l) {
                    pop_l -= (((*cigar)[i]) >> 4);
                    (*cigar_n) -= 1;
                } else if ((((*cigar)[i]) >> 4) == pop_l) {
                    (*cigar_n) -= 1;
                    return 0;
                } else {
                    long long op = (*cigar)[i] & 0xf;
                    (*cigar)[i] -= ((pop_l << 4) | op);
                    return 0;
                }
            } 
            //else if ((((*cigar)[i]) & 0xf) == CINS)
            //else if ((((*cigar)[i]) & 0xf) == CDEL)
            --i;
            if (i < 0) { fprintf(stderr, "\n[pop_cigar] Cigar error,\n"); exit(1); }
        }
    } else {
        fprintf(stderr, "\n[pop_cigar] Unknown pop mode. (%d)\n", pop_mode); exit(1);
    }
    return 0;
}

static inline int _invert_cigar(cigar32_t **cigar, int cigar_n) {
    if (cigar_n <= 1) return 0;
    int i;
    cigar32_t tmp;
    for (i = 0; i < cigar_n/2; ++i) {
        tmp = (*cigar)[i];
        (*cigar)[i] = (*cigar)[cigar_n - i - 1];
        (*cigar)[cigar_n - i - 1] = tmp;
    }
    return 0;
}

static inline int _push_cigar0(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t _cigar) {
    int i;
    i = *cigar_n;
    if (((i-1) >=0) && (((*cigar)[i-1] & 0xf) == (_cigar & 0xf)))
        (*cigar)[i-1] += ((_cigar >> 4) << 4);
    else {
        if (i == *cigar_m) {
            (*cigar_m) = (*cigar_m) ? (*cigar_m)<<1 : 4;
			(*cigar) = (cigar32_t*)realloc(*cigar, (*cigar_m) * sizeof (cigar32_t));
			if ((*cigar) == NULL)	{fprintf(stderr, "\n[frag_check] Memory is not enougy.\n");exit(1);}
        }
        (*cigar)[i] = _cigar;
        ++(*cigar_n);
    }
    return 0;
}

static inline int _push_cigar1(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t _cigar) {
	if ((_cigar >> 4) == 0) return 0;
	return _push_cigar0(cigar, cigar_n, cigar_m, _cigar);
}

static inline int _push_cigar(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t *_cigar, int _cigar_n) {
    if (_cigar_n == 0) return 0;
	int i,j;
	i = *cigar_n;

    if ((i-1) >= 0) {
        if(((*cigar)[i-1] & 0xf) == (_cigar[0] & 0xf)) {
            (*cigar)[i-1] += ((_cigar[0] >> 4) << 4);
            j = 1;
        } else if ((((*cigar)[i-1] & 0xf) == CINS && (_cigar[0] & 0xf) == CSOFT_CLIP)   // xI+yS -> (x+y)S 
                ||(((*cigar)[i-1] & 0xf) == CSOFT_CLIP && (_cigar[0] & 0xf) == CINS)){ // xS+yI -> (x+y)S 
            (*cigar)[i-1] = ((((*cigar)[i-1]>>4) + (_cigar[0]>>4)) << 4) | CSOFT_CLIP;
            j = 1;
        } else j = 0;
    } else j = 0;

    for (; j < _cigar_n; ++i,++j) {
        if (i == *cigar_m) {
            (*cigar_m) = (*cigar_m) ? (*cigar_m)<<1 : 4;
            (*cigar) = (cigar32_t*)realloc(*cigar, (*cigar_m) * sizeof (cigar32_t));
			if ((*cigar) == NULL)	{fprintf(stderr, "\n[frag_check] Memory is not enougy.\n");exit(1);}
		}
		(*cigar)[i] = _cigar[j];
	}
	*cigar_n = i;
    return 0;
}

static inline int _push_cigar_e1(cigar32_t **cigar, int *cigar_n, int *cigar_m, uint64_t *refend, int *readend, cigar32_t _cigar) {
    _push_cigar1(cigar, cigar_n, cigar_m, _cigar);
    if (refend) *refend += refInCigar(&_cigar, 1);
    if (readend) *readend += readInCigar(&_cigar, 1);
    return 0;
}

static inline int _push_cigar_e(cigar32_t **cigar, int *cigar_n, int *cigar_m, uint64_t *refend, int *readend, cigar32_t *_cigar, int _cigar_n) {
    _push_cigar(cigar, cigar_n, cigar_m, _cigar, _cigar_n);
    if (refend) *refend += refInCigar(_cigar, _cigar_n);
    if (readend) *readend += readInCigar(_cigar, _cigar_n);
    return 0;
}


//int _push_cigar(cigar32_t **cigar, int *cigar_len, int *cigar_m, cigar32_t *_cigar, int _cigar_len);
//int _push_cigar1(cigar32_t **cigar, int *cigar_len, int *cigar_m, cigar32_t _cigar);
//int _push_cigar0(cigar32_t **cigar, int *cigar_len, int *cigar_m, cigar32_t _cigar);
//void _invert_cigar(cigar32_t **cigar, int cigar_n);
void frag_init_msg(frag_msg *f_msg, int frag_max);
void frag_free_msg(frag_msg *f_msg, int line_num);
int frag_set_msg(map_msg *m_msg, int seed_i, int aln_i, int FLAG, frag_msg *f_msg, int frag_i);//FLAG 0: start/1:end / 2:seed
int frag_copy_msg(frag_msg *ff_msg, frag_msg *tf_msg);
void lamsa_res_aux(line_aln_res *la, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, int read_len, lamsa_aln_para *AP, kseq_t *seqs);
void frag_check(map_msg *m_msg, frag_msg **f_msg, aln_res *a_res,
               bntseq_t *bns, uint8_t *pac, 
               uint8_t *read_bseq, uint8_t **read_rbseq,
               lamsa_aln_per_para *APP, lamsa_aln_para *AP,
               kseq_t *seqs,
               int line_n,
               uint32_t **hash_num, uint64_t ***hash_node);

void check_cigar(cigar32_t *cigar, int cigar_len, char *read_name, int read_len);
void printcigar(FILE *outp, cigar32_t *cigar, int cigar_len);

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
