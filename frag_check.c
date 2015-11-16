/*
 * frag_check.c
 *
 * function: transform raw frag.msg to detailed frag.msg
 * method:   extend-ssw
 *
 * Created by Yan Gao on 08/29/2013
 */
/*
 * Convert the whole read seq to its rev-com seq when its aln-result is '-' srand.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <time.h>
#include "frag_check.h"
#include "lsat_aln.h"
#include "split_mapping.h"
#include "bntseq.h"
#include "kseq.h"
#include "ksw.h"

#define PER_LEN 100

KSEQ_INIT(gzFile, gzread)
const int8_t bwasw_sc_mat[25] = {1, -3, -3, -3, -1,
							    -3,  1, -3, -3, -1,
							    -3, -3,  1, -3, -1,
							    -3, -3, -3,  1, -1,
							    -1, -1, -1, -1,  1};

const int8_t sc_mat[25] =       {1, -2, -2, -2, -1,
				                -2,  1, -2, -2, -1,
				                -2, -2,  1, -2, -1,
				                -2, -2, -2,  1, -1,
				                -1, -1, -1, -1, -1};

const int8_t ssw_sc_mat[25] =   {2, -2, -2, -2, -1,
							    -2,  2, -2, -2, -1,
							    -2, -2,  2, -2, -1,
							    -2, -2, -2,  2, -1,
							    -1, -1, -1, -1,  1};
/*
 * parameter: frag.msg
 */
// for debug
char READ_NAME[1024];

void frag_init_msg(frag_msg *f_msg, int frag_max)
{
	f_msg->frag_max = frag_max; 
	f_msg->frag_num = 0;
	f_msg->per_seed_max = frag_max;
	//f_msg->seed_num = 0;

	f_msg->fa_msg = (frag_aln_msg*)malloc(frag_max*sizeof(frag_aln_msg));
	int i;
	for (i = 0; i < frag_max; ++i)
	{
		f_msg->fa_msg[i].cigar_max = 100;
		f_msg->fa_msg[i].cigar = (cigar32_t*)malloc(100 * sizeof(cigar32_t));
		f_msg->fa_msg[i].cigar_len = 0;

		f_msg->fa_msg[i].seed_i = (int*)malloc(f_msg->per_seed_max*sizeof(int));
		f_msg->fa_msg[i].seed_aln_i = (int*)malloc(f_msg->per_seed_max*sizeof(int));
		f_msg->fa_msg[i].seed_num = 0;
	}
}

int frag_copy_msg(frag_msg *ff_msg, frag_msg *tf_msg)
{
	tf_msg->frag_max = ff_msg->frag_max;
	tf_msg->frag_num = 0;
	//tf_msg->last_len = ff_msg->last_len;
	//tf_msg->seed_num = 0;
	//tf_msg->seed_all = ff_msg->seed_all;
	tf_msg->per_seed_max = ff_msg->per_seed_max;
	//tf_msg->fa_msg = (frag_aln_msg*)malloc(tf_msg->frag_max * sizeof(frag_aln_msg));
	int i;
	for (i = 0; i < tf_msg->frag_max; ++i) {
		tf_msg->fa_msg[i].cigar_max = 100;
		//tf_msg->fa_msg[i].cigar = (ucigar32_t*)malloc(100 * sizeof(ucigar32_t));
		tf_msg->fa_msg[i].cigar_len = 0;

		//tf_msg->fa_msg[i].seed_i = (int*)malloc(tf_msg->per_seed_max*sizeof(int));
		//tf_msg->fa_msg[i].seed_aln_i = (int*)malloc(tf_msg->per_seed_max*sizeof(int));
		tf_msg->fa_msg[i].seed_num = 0;
	}
	return 0;
}

void frag_free_msg(frag_msg *f_msg, int line_num)
{
	int i, j;
	for (i = 0; i < line_num; ++i)
	{
        frag_msg *p = f_msg+i;
		for (j = 0; j < p->frag_max; ++j) {
			free(p->fa_msg[j].cigar);
			free(p->fa_msg[j].seed_i);
			free(p->fa_msg[j].seed_aln_i);
		}
		free(p->fa_msg);
	}
	free(f_msg);
}

int frag_set_msg(aln_msg *a_msg, int seed_i, int aln_i,
				 int FLAG, frag_msg *f_msg, int frag_i)//FLAG 0:start / 1:end / 2:seed
{
	if (FLAG == FRAG_END) {	//end
		f_msg->fa_msg[frag_i].chr = a_msg[seed_i].at[aln_i].chr;
		f_msg->fa_msg[frag_i].srand = a_msg[seed_i].at[aln_i].nsrand;	//+:1/-:-1
		f_msg->fa_msg[frag_i].seed_i[0] = seed_i;
		f_msg->fa_msg[frag_i].seed_aln_i[0] = aln_i;
		f_msg->fa_msg[frag_i].seed_num = 1;
		f_msg->fa_msg[frag_i].flag = UNCOVERED;
		f_msg->fa_msg[frag_i].cigar_len = 0;
	}
	else if(FLAG==FRAG_START) {	//start
		f_msg->frag_num = frag_i + 1;
	}
	else {			//seed
		f_msg->fa_msg[frag_i].seed_i[f_msg->fa_msg[frag_i].seed_num] = seed_i;
		f_msg->fa_msg[frag_i].seed_aln_i[f_msg->fa_msg[frag_i].seed_num] = aln_i;
		f_msg->fa_msg[frag_i].seed_num++;
	}
	return 0;
}

int get_ref_intv(uint8_t **ref_bseq, int *ref_max_blen, 
				 bntseq_t *bns, uint8_t *pac, 
				 aln_msg *a_msg, 
				 int seed1_i, int seed1_aln_i, int seed2_i, int seed2_aln_i, 
                 lsat_aln_per_para APP)
{
	int64_t start; int32_t len;
    start = a_msg[seed1_i].at[seed1_aln_i].offset + APP.seed_len - 1 + a_msg[seed1_i].at[seed1_aln_i].len_dif;
    len = a_msg[seed2_i].at[seed2_aln_i].offset - 1 - start;
    if (len <= 0) return 0;
	if (len > *ref_max_blen) {
		*ref_bseq = (uint8_t*)realloc(*ref_bseq, len * sizeof(uint8_t));
		*ref_max_blen = len;
	}
    //printf("offset: %lld ", (long long)start);
    pac2fa_core(bns, pac, a_msg[seed1_i].at[seed1_aln_i].chr, start, &len, *ref_bseq);
    //int j;
	//fprintf(stdout, "ref:\t%d\t%ld\t", len, start);
	//for (j = 0; j < len; j++)
	//	fprintf(stdout, "%d", (int)seq1[j]);
	//fprintf(stdout, "\n");
	if (len > 100)
		fprintf(stderr, "debug");
	return (int)len;
}

int get_read_intv(uint8_t *seq2, uint8_t *read_bseq, 
				  aln_msg *a_msg, 
				  int seed1_i, int seed1_aln_i, int seed2_i, int seed2_aln_i, 
				  int *band_width, 
				  lsat_aln_per_para APP)
{
	int32_t i;
	int j;

	//set band-width
	*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
	
	
	//convert char seq to int seq

	if (a_msg[seed1_i].at[seed1_aln_i].nsrand == 1)
	{
		*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
		for (j=0, i = a_msg[seed1_i].read_id * APP.seed_step - APP.seed_inv; i < (a_msg[seed2_i].read_id-1) * APP.seed_step; ++j, ++i)
			seq2[j] = read_bseq[i];
	}
	else
	{
		//printf("read: %d, %d ", a_msg[seed1_i].read_id,a_msg[seed1_i].read_id * 2 * seed_len);
		*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
        for (j=0, i = APP.last_len+a_msg[seed1_i].read_id*APP.seed_step-APP.seed_inv; i < APP.last_len+(a_msg[seed2_i].read_id-1)*APP.seed_step; ++j, ++i)
		//for (j=0, i = a_msg[seed1_i].read_id * APP.seed_step; i < (a_msg[seed2_i].read_id * 2 - 1) * seed_len; ++j, ++i)
		{
			seq2[j] = read_bseq[i];
			//printf("%c", read_seq[i-100]);
		}
	}
	/*else 
	{
		*band_width = 2 * (((a_msg[seed1_i].read_id - a_msg[seed2_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
		read_len = (((a_msg[seed1_i].read_id - a_msg[seed2_i].read_id) << 1) - 1) * seed_len;
		for (j=0; j < read_len; ++j)
			seq2[read_len-j-1] = com_nst_nt4_table[(int)read_seq[((a_msg[seed2_i].read_id << 1) - 1) * seed_len + j]];
	}*/
	
	//fprintf(stdout, "read:\t%d\t%d\t", j, s);
	//for (i = 0; i < j; i++)
	//	fprintf(stdout, "%d",(int)seq2[i]);
	//fprintf(stdout, "\n");
	return (j);
}

void printcigar(FILE *outp, cigar32_t *cigar, int cigar_len)
{
	int i;
	for (i = 0; i < cigar_len; i++)
		fprintf(outp, "%d%c", (int)(cigar[i]>>4), CIGAR_STR[(int)(cigar[i] & 0xf)]);
}

void printnst(char *msg1, uint8_t *seq, int len, char *msg2)
{
	int i;
	printf("%s", msg1);
	for (i = 0; i < len; ++i) printf("%c", "ACGTN"[(int)seq[i]]);
	printf("%s", msg2);
}

cigar_t* _init_cigar(int cigar_n)
{
    cigar_t *ct = (cigar_t*)malloc(sizeof(cigar_t));
    ct->cigar_n = cigar_n;
    ct->cigar_m = cigar_n;
    ct->cigar = (cigar32_t*)malloc(sizeof(cigar32_t));
    return ct;
}
void _free_cigar(cigar_t *ct)
{
    free(ct->cigar);
    free(ct);
}

//return 'MIS' length in cigar
int readInCigar(cigar32_t *cigar, int cigar_len) {
    int i, read_len=0;
    for (i = 0; i < cigar_len; ++i) {
        if ((cigar[i] & 0xf) == CMATCH || (cigar[i] & 0xf) == CINS || (cigar[i] & 0xf) == CSOFT_CLIP)
            read_len += (cigar[i]  >> 4);
        else if ((cigar[i] & 0xf) != CDEL && (cigar[i] & 0xf) != CHARD_CLIP) {
            fprintf(stderr, "\n%s\n[readInCigar] Cigar Error.\n", READ_NAME); printcigar(stderr, cigar, cigar_len); exit(1);
        }
    }
    return read_len;
}
//
//return 'MDH' length in cigar
int refInCigar(cigar32_t *cigar, int cigar_len) {
    int i, read_len=0;
    for (i = 0; i < cigar_len; ++i) {
        if ((cigar[i] & 0xf) == CMATCH || (cigar[i] & 0xf) == CDEL || (cigar[i] & 0xf) == CHARD_CLIP)
            read_len += (cigar[i]  >> 4);
        else if ((cigar[i] & 0xf) != CINS && (cigar[i] & 0xf) != CSOFT_CLIP) {
            fprintf(stderr, "\n%s\n[readInCigar] Cigar Error.\n", READ_NAME); printcigar(stderr, cigar, cigar_len); exit(1);
        }
    }
    return read_len;
}

void _deQueue_cigar(cigar32_t **cigar, int *cigar_n, int de_n) {
    if (de_n == 0) return;
    if (de_n >= *cigar_n) { *cigar_n = 0; return; }

    int i;
    for (i = de_n; i < *cigar_n; ++i)
        (*cigar)[i - de_n] = (*cigar)[i];
    (*cigar_n) -= de_n;
}

//pop_mode: 0 pop_n -> cigar_n
//          1 pop_n -> base n
void _pop_cigar(cigar32_t **cigar, int *cigar_n, int pop_mode, int pop_n) {
    if (pop_mode == 0)
    {
        (*cigar_n) -= pop_n;
        if (*cigar_n < 0) *cigar_n = 0;
    }
    else if (pop_mode == 1)
    {
        int i;
        int pop_l = pop_n;

        i = *cigar_n - 1;
        while (pop_l != 0)
        {
            if ((((*cigar)[i]) & 0xf) == CMATCH || (((*cigar)[i]) & 0xf) == CINS)
            {
                if ((((*cigar)[i]) >> 4) < pop_l)
                {
                    pop_l -= (((*cigar)[i]) >> 4);
                    (*cigar_n) -= 1;
                }
                else if ((((*cigar)[i]) >> 4) == pop_l)
                {
                    (*cigar_n) -= 1;
                    return;
                }
                else
                {
                    long long op = (*cigar)[i] & 0xf;
                    (*cigar)[i] -= ((pop_l << 4) | op);
                    return;
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
}

void _invert_cigar(cigar32_t **cigar, int cigar_n) {
    if (cigar_n <= 1) return;
    int i;
    cigar32_t tmp;
    for (i = 0; i < cigar_n/2; ++i) {
        tmp = (*cigar)[i];
        (*cigar)[i] = (*cigar)[cigar_n - i - 1];
        (*cigar)[cigar_n - i - 1] = tmp;
    }
}

void _push_cigar0(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t _cigar) {
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
}

void _push_cigar1(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t _cigar) {
	if ((_cigar >> 4) == 0) return;
	_push_cigar0(cigar, cigar_n, cigar_m, _cigar);
}

void _push_cigar(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t *_cigar, int _cigar_n)
{
    if (_cigar_n == 0) return;
	int i,j;
	i = *cigar_n;
	if (((i-1) >= 0) && (((*cigar)[i-1] & 0xf) == (_cigar[0] & 0xf))) {
		(*cigar)[i-1] += ((_cigar[0] >> 4) << 4);
		j = 1;
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
}

void _push_cigar_e1(cigar32_t **cigar, int *cigar_n, int *cigar_m, uint64_t *refend, int *readend, cigar32_t _cigar) {
    _push_cigar1(cigar, cigar_n, cigar_m, _cigar);
    if (refend) *refend += refInCigar(&_cigar, 1);
    if (readend) *readend += readInCigar(&_cigar, 1);
}

void _push_cigar_e(cigar32_t **cigar, int *cigar_n, int *cigar_m, uint64_t *refend, int *readend, cigar32_t *_cigar, int _cigar_n) {
    _push_cigar(cigar, cigar_n, cigar_m, _cigar, _cigar_n);
    if (refend) *refend += refInCigar(_cigar, _cigar_n);
    if (readend) *readend += readInCigar(_cigar, _cigar_n);
}

//chr/nsrand
void push_res(line_aln_res *la)
{
    int i;
    if (la->cur_res_n == la->res_m-1) {
        la->res_m <<= 1;
        la->res = (res_t*)realloc(la->res, la->res_m * sizeof(res_t));
        if (la->res == NULL) { fprintf(stderr, "[lsat_res_split] Not enough memory.\n"); exit(0); }
        for (i = la->cur_res_n+1; i < la->res_m; ++i) {
            la->res[i].c_m = 10;
            la->res[i].cigar = (cigar32_t*)malloc(10 * sizeof(cigar32_t));
            la->res[i].cigar_len = 0;
        }
    }
    ++(la->cur_res_n);
    la->res[la->cur_res_n].chr = la->res[la->cur_res_n-1].chr;
    la->res[la->cur_res_n].nsrand = la->res[la->cur_res_n-1].nsrand;
}

//merge interval and seed to frag
//frag has the first seed's msg already
void merge_cigar(cigar32_t **c1, int *c1_n, int *c1_m, uint64_t *c1_refend, int *c1_readend, int chr,
                  cigar32_t *_c2, int c2_n, int c2_reflen, int c2_readlen, 
				  bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, lsat_aln_para AP)
{
    if (c2_n == 0) return;
	if (*c1_n !=0  && (((((*c1)[*c1_n-1] & 0xf) == CINS || ((*c1)[*c1_n-1] & 0xf) == CDEL) && (_c2[0] & 0xf) != CSOFT_CLIP && (_c2[0] & 0xf) != CHARD_CLIP)
	   || (((_c2[0] & 0xf) == CINS || (_c2[0] & 0xf) == CDEL) && ((*c1)[*c1_n-1] & 0xf) != CSOFT_CLIP && ((*c1)[*c1_n-1] & 0xf) != CHARD_CLIP))) { // repair boundary
		    /* seq1, ref */ /*    seq2, read     */
		int len1, len11=0,  len2, len21=0, len22=0, len_dif1=0, len_dif2=0;
		int bd_cigar_len=0, b=0, min_b, ci=0; cigar32_t *bd_cigar = 0;
		int i, j, left=1, right=1;
		cigar32_t *c2 = (cigar32_t*)malloc(c2_n * sizeof(cigar32_t));
		for (i = 0; i < c2_n; ++i) c2[i] = _c2[i];
		uint8_t *seq1, *seq2;
        int64_t ref_start; int read_start;
		int md = 5; 
		while (1) {	//get a right boundary-cigar(without 'I'/'D' on the head or tail)
			//extend by md-bp
			if (left) {
				while ((*c1_n) >= 1) {	// get a nM, n > 3
					if (((*c1)[*c1_n-1] & 0xf) == CMATCH && ((*c1)[*c1_n-1]>>4) > md) {
						(*c1)[*c1_n-1] -= (md << 4); len21 += md; break;}
					else if (((*c1)[*c1_n-1] & 0xf) == CMATCH) {
						len21 += (int)((*c1)[*c1_n-1] >> 4); --(*c1_n); }
					else if (((*c1)[*c1_n-1] & 0xf) == CINS) {
						len21 += (int)((*c1)[*c1_n-1] >> 4); len_dif1 -= (int)((*c1)[*c1_n-1] >> 4); b += (int)((*c1)[*c1_n-1]>>4); --(*c1_n); }
					else if (((*c1)[*c1_n-1] & 0xf) == CDEL) {
						len_dif1 += (int)((*c1)[*c1_n-1] >> 4); b += (int)((*c1)[*c1_n-1] >> 4); --(*c1_n); }
					else { 
						left = -1; break;
						fprintf(stderr, "\n%s\n[merge cigar] Unexpected(1): ", READ_NAME); 
						printcigar(stderr, (*c1), (*c1_n)); exit(1);
					}
				}
				len11 = len21 + len_dif1;
				read_start = *c1_readend - len21 + 1; // 1-based
				ref_start = *c1_refend - len11 + 1;	    
			}
			//extend 3-match boundary
			if (right) {
				while (ci < c2_n) {	// get a nM, n > 3
					if ((c2[ci] & 0xf) == CMATCH && (c2[ci] >> 4) > md) {
						c2[ci] -= (md << 4); len22 += md; break;}
					else if ((c2[ci] & 0xf) == CMATCH) {
						len22 += (int)(c2[ci] >> 4); ci++; }
					else if ((c2[ci] & 0xf) == CINS) {
						len22 += (int)(c2[ci] >> 4); len_dif2 -= (int)(c2[ci]>>4); b += (int)(c2[ci]>>4); ++ci; }
					else if ((c2[ci] & 0xf) == CDEL) {
						len_dif2 += (int)(c2[ci] >> 4); b += (int)(c2[ci]>>4); ++ci; }
					else { 
						right = -1; break;
						fprintf(stderr, "\n%s\n[merge cigar] Unexpected(2): ", READ_NAME); 
						printcigar(stderr, c2+ci, c2_n-ci); exit(1);
					}
				}
			}
			len2 = len21+len22; len1 = len2 + len_dif1+len_dif2;
			min_b = abs(len_dif1+len_dif2) + md; b = MAXOFTWO(b, min_b);
			seq1 = (uint8_t *)malloc((len1+1) * sizeof(uint8_t)); seq2 = (uint8_t *)malloc((len2+1) * sizeof(uint8_t));
			pac2fa_core(bns, pac, chr, ref_start-1, &len1, seq1);
			for (i=0, j= read_start-1; j<read_start+len2-1; ++i, ++j) {
				seq2[i] = read_bseq[j];
			}
			ksw_global(len2, seq2, len1, seq1, 5, bwasw_sc_mat, AP.gapo, AP.gape, b, &bd_cigar_len, &bd_cigar);
			free(seq1); free(seq2);

			if ((bd_cigar[0] & 0xf) == CMATCH) left = 0;
			else if (*c1_n == 0 || left < 0) break;  // unfinished merge
			if ((bd_cigar[bd_cigar_len-1] & 0xf) == CMATCH) right = 0;
			else if (ci == c2_n || right < 0) break; // unfinished merge

			if ((left + right) == 0 ) break;  // finished merge
			//left = right = 1;
			free(bd_cigar); 
		}
		_push_cigar(c1, c1_n, c1_m, bd_cigar, bd_cigar_len);
		_push_cigar(c1, c1_n, c1_m, c2+ci, c2_n-ci);
		free(bd_cigar); if (c2) free(c2);
	} else _push_cigar(c1, c1_n, c1_m, _c2, c2_n);
    (*c1_refend) +=  c2_reflen;
    (*c1_readend) +=  c2_readlen;
}

//Extend the intervals between seeds in the same frag
//un-seed region: end of + srand and start of - end 
int frag_extend(frag_msg *f_msg, aln_msg *a_msg, int f_i, 
				bntseq_t *bns, uint8_t *pac, 
				uint8_t *read_bseq, uint8_t *bseq1, uint8_t **bseq2, int *ref_max_blen,
                lsat_aln_per_para APP, lsat_aln_para AP,
				line_aln_res *la)
{
    frag_aln_msg fa_msg = f_msg->fa_msg[f_i];
	int i, seed_i, seed_aln_i, last_i, last_aln_i;
	int len1, len2; /*read, ref*/
	int b, min_b, cigar_len;

	//banded global alignment for intervals between seeds
	cigar_len = 0;
	
	//get the first node
	if (fa_msg.srand == 1) {
		i = fa_msg.seed_num-1;
		last_i = fa_msg.seed_i[i];
		last_aln_i = fa_msg.seed_aln_i[i];
        fa_msg.cigar_read_start = (a_msg[last_i].read_id-1) * APP.seed_step + 1;
        fa_msg.cigar_read_end = fa_msg.cigar_read_start - 1 + APP.seed_len;
	} else {	//'-' srand 
		last_i = fa_msg.seed_i[0];
		last_aln_i = fa_msg.seed_aln_i[0];
        fa_msg.cigar_read_start = APP.last_len + (a_msg[last_i].read_id-1) * APP.seed_step + 1; 
        fa_msg.cigar_read_end = fa_msg.cigar_read_start - 1 + APP.seed_len;
	}

	//copy the first seed's cigar to frag
	{
		_push_cigar(&(fa_msg.cigar), &(fa_msg.cigar_len), &(fa_msg.cigar_max), a_msg[last_i].at[last_aln_i].cigar, a_msg[last_i].at[last_aln_i].cigar_len);
		//start/end: 1-base
		fa_msg.cigar_ref_start = a_msg[last_i].at[last_aln_i].offset;
		fa_msg.cigar_ref_end = a_msg[last_i].at[last_aln_i].offset + APP.seed_len -1 + a_msg[last_i].at[last_aln_i].len_dif;
	}
	if (fa_msg.srand == 1)
	{
		for (--i; i >= 0; --i) {
			seed_i = fa_msg.seed_i[i];
			seed_aln_i = fa_msg.seed_aln_i[i];
			//get ref seq and seed-interval seq
			if ((len2 = get_ref_intv(bseq2, ref_max_blen, bns, pac, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, APP)) < 0) { fprintf(stderr, "\n[frag extend] frag path error : (%d,%d) -> %d %lld %d %lld\n", last_i, last_aln_i, a_msg[last_i].read_id, (long long)a_msg[last_i].at[last_aln_i].offset, a_msg[seed_i].read_id, (long long)a_msg[seed_i].at[seed_aln_i].offset); return -1;}
			len1 = get_read_intv(bseq1, read_bseq, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, &b, APP);
			min_b = abs(len2-len1)+3;
			b = MAXOFTWO(b, min_b);
			cigar_len = 0;
			cigar32_t *cigar=0;
			ksw_global(len1, bseq1, len2, *bseq2, 5, bwasw_sc_mat, AP.gapo, AP.gape, b, &cigar_len, &cigar);
			merge_cigar(&(fa_msg.cigar), &(fa_msg.cigar_len), &(fa_msg.cigar_max), &(fa_msg.cigar_ref_end), &(fa_msg.cigar_read_end), fa_msg.chr,
					     cigar, cigar_len, len2, len1, bns, pac, read_bseq, AP);
			merge_cigar(&(fa_msg.cigar), &(fa_msg.cigar_len), &(fa_msg.cigar_max), &(fa_msg.cigar_ref_end), &(fa_msg.cigar_read_end), fa_msg.chr, 
					     a_msg[seed_i].at[seed_aln_i].cigar, a_msg[seed_i].at[seed_aln_i].cigar_len, APP.seed_len+a_msg[seed_i].at[seed_aln_i].len_dif, APP.seed_len, bns, pac, read_bseq, AP);//merge seed to frag

			last_i = seed_i;
			last_aln_i = seed_aln_i;
			free(cigar);
		}
	} else {
		for (i = 1; i < fa_msg.seed_num; ++i)
		{
			seed_i = fa_msg.seed_i[i];
			seed_aln_i = fa_msg.seed_aln_i[i];
			//get ref seq and seed-interval seq
			if ((len2 = get_ref_intv(bseq2, ref_max_blen, bns, pac, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, APP)) < 0) { fprintf(stderr, "\n[frag extend] frag path error : (%d,%d) -> %d %lld %d %lld\n", last_i, last_aln_i, a_msg[last_i].read_id, (long long)a_msg[last_i].at[last_aln_i].offset, a_msg[seed_i].read_id, (long long)a_msg[seed_i].at[seed_aln_i].offset); return -1;}
			len1 = get_read_intv(bseq1, read_bseq, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, &b, APP);
			min_b = abs(len2-len1)+3;
			b = MAXOFTWO(b, min_b);
			cigar_len = 0;
			cigar32_t *cigar=0;
			ksw_global(len1, bseq1, len2, *bseq2, 5, bwasw_sc_mat, AP.gapo, AP.gape, b, &cigar_len, &cigar);
            merge_cigar(&(fa_msg.cigar), &(fa_msg.cigar_len), &(fa_msg.cigar_max), &(fa_msg.cigar_ref_end), &(fa_msg.cigar_read_end), fa_msg.chr,
                         cigar, cigar_len, len2, len1, bns, pac, read_bseq, AP);
            merge_cigar(&(fa_msg.cigar), &(fa_msg.cigar_len), &(fa_msg.cigar_max), &(fa_msg.cigar_ref_end), &(fa_msg.cigar_read_end), fa_msg.chr,
                        a_msg[seed_i].at[seed_aln_i].cigar, a_msg[seed_i].at[seed_aln_i].cigar_len, APP.seed_len+a_msg[seed_i].at[seed_aln_i].len_dif, APP.seed_len, bns, pac, read_bseq, AP);//merge seed to frag 

			last_i = seed_i;
			last_aln_i = seed_aln_i;
			free(cigar);
		}
	}
    merge_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), &(la->res[la->cur_res_n].refend), &(la->res[la->cur_res_n].readend), fa_msg.chr, fa_msg.cigar, fa_msg.cigar_len, fa_msg.cigar_ref_end-fa_msg.cigar_ref_start+1, fa_msg.cigar_read_end-fa_msg.cigar_read_start+1, bns, pac, read_bseq, AP);
	return 0;
}

//for '-' srand:
//1. convert read seq into reverse-complement seq.
//2. split-map the re-co seq to the ref.
//3. combine with aln results of other regions.
void split_mapping(bntseq_t *bns, uint8_t *pac, 
				   uint8_t *read_bseq, 
				   frag_msg *f_msg, aln_msg *a_msg, 
				   uint32_t **hash_num, uint64_t ***hash_node, 
				   int f1_i, int f2_i, 
                   lsat_aln_per_para APP, lsat_aln_para AP,
				   line_aln_res *la)
{
	cigar32_t *s_cigar = (cigar32_t*)malloc(2 * sizeof(cigar32_t));
    int s_clen = 0, s_cm = 2;
    cigar32_t *_cigar=0; int _cigar_n=0, _cigar_m;
    int i;
    int nsrand, s1_i, s1_aln_i, seed_num, s2_i, s2_aln_i;
    frag_aln_msg fa_msg1 = f_msg->fa_msg[f1_i];
    frag_aln_msg fa_msg2 = f_msg->fa_msg[f2_i];

	nsrand = fa_msg1.srand;
	if (nsrand == 1) {
		s1_i = fa_msg1.seed_i[0];
		s1_aln_i = fa_msg1.seed_aln_i[0];

		seed_num = fa_msg2.seed_num;
		s2_i = fa_msg2.seed_i[seed_num-1];
		s2_aln_i = fa_msg2.seed_aln_i[seed_num-1];
	} else {
		seed_num = fa_msg1.seed_num;
		s1_i = fa_msg1.seed_i[seed_num-1];
		s1_aln_i = fa_msg1.seed_aln_i[seed_num-1];

		s2_i = fa_msg2.seed_i[0];
		s2_aln_i = fa_msg2.seed_aln_i[0];
	}
    aln_t at1 = a_msg[s1_i].at[s1_aln_i];
    aln_t at2 = a_msg[s2_i].at[s2_aln_i];

	int s_qlen, s_tlen, res;
	uint8_t *s_qseq, *s_tseq;
	int64_t ref_offset;
	//init value for hash-map
	int hash_len = AP.hash_len;

	s_qlen = (a_msg[s2_i].read_id-a_msg[s1_i].read_id) * APP.seed_step - APP.seed_len;
	s_qseq = (uint8_t*)malloc(s_qlen * sizeof(uint8_t));

    int bd;
    get_read_intv(s_qseq, read_bseq, a_msg, s1_i, s1_aln_i, s2_i, s2_aln_i, &bd, APP);

	//check SV-type
	{
		int64_t pos = at1.offset+APP.seed_len-1+at1.len_dif;
		int64_t exp = at1.offset + at1.len_dif + (a_msg[s2_i].read_id - a_msg[s1_i].read_id) * APP.seed_step;	
		int64_t act = at2.offset;
		int dis = act-exp;

		//fprintf(stdout, "s1: %lld, exp_s2: %lld act_s2: %lld %d\n", (long long)pos, (long long)exp, (long long)act, dis);
#ifdef __DEBUG__
        if (abs(dis) > APP.match_dis) 
			fprintf(stdout, "%d\t%lld\t%d\t%s\n", at1.chr, (long long)pos, abs(dis), dis>0?"DEL":"INS");
#endif
		if (dis > APP.match_dis) {	//DEL
            s_tlen = s_qlen + dis;
            s_tseq = (uint8_t*)malloc(s_tlen * sizeof(uint8_t));
            ref_offset = at1.offset + APP.seed_len + at1.len_dif;
            pac2fa_core(bns, pac, at1.chr, ref_offset-1, &s_tlen, s_tseq);
            if (s_qlen < hash_len) {
                ksw_bi_extend(s_qlen, s_qseq, s_tlen, s_tseq, 5, bwasw_sc_mat, dis+3, hash_len*bwasw_sc_mat[0], hash_len*bwasw_sc_mat[0], AP, &_cigar, &_cigar_n, &_cigar_m);
                _push_cigar(&s_cigar, &s_clen, &s_cm, _cigar, _cigar_n);
                free(_cigar);
            } else {
                res = split_indel_map(&s_cigar, &s_clen, &s_cm, s_qseq, s_qlen, s_tseq, s_tlen, 0, AP,  hash_num, hash_node);
            }
        } else if (dis < -APP.match_dis) { 
            s_tlen = s_qlen + dis;
            if (s_tlen < 0) { // overlapped ins
                int _s_tlen = s_qlen + hash_len;
                s_tseq = (uint8_t*)malloc(_s_tlen * sizeof(uint8_t));
                int lqe, lte, rqe, rte;

                // left-extend
                ref_offset = at1.offset + APP.seed_len + at1.len_dif;
                pac2fa_core(bns, pac, at1.chr, ref_offset-1, &_s_tlen, s_tseq);
                ksw_extend_core(s_qlen, s_qseq, _s_tlen, s_tseq, 5, bwasw_sc_mat, 3, hash_len*bwasw_sc_mat[0], AP, &lqe, &lte, &_cigar, &_cigar_n, &_cigar_m);
                _push_cigar(&s_cigar, &s_clen, &s_cm, _cigar, _cigar_n);
                free(_cigar);

                // right-extend
                ref_offset = at2.offset - _s_tlen;
                pac2fa_core(bns, pac, at2.chr, ref_offset-1, &_s_tlen, s_tseq);
                // reverse ref and read seq
				char tmp;
                for (i = 0; i < s_qlen>>1; ++i) {
					tmp = s_qseq[s_qlen-1-i];
					s_qseq[s_qlen-1-i] = s_qseq[i];
					s_qseq[i] = tmp;
				}
                for (i = 0; i < _s_tlen>>1; ++i) {
					tmp = s_tseq[_s_tlen-1-i];
					s_tseq[_s_tlen-1-i] = s_tseq[i];
					s_tseq[i] = tmp;
				}
                ksw_extend_core(s_qlen, s_qseq, _s_tlen, s_tseq, 5, bwasw_sc_mat, 3, hash_len*bwasw_sc_mat[0], AP, &rqe, &rte, &_cigar, &_cigar_n, &_cigar_m);
                _invert_cigar(&_cigar, _cigar_n);
                
                // merge, add overlap-flag('S/H')
                extern void sw_mid_fix(cigar32_t **cigar, int *cigar_n, int *cigar_m, const uint8_t *query, int qlen, int lqe, int rqe, cigar32_t *rcigar, int rn_cigar, const uint8_t *target, int tlen, int lte, int rte, lsat_aln_para AP, int m, const int8_t *mat);
                sw_mid_fix(&s_cigar, &s_clen, &s_cm, s_tseq, s_qlen, lqe, rqe, _cigar, _cigar_n, s_tseq, s_qlen+dis, lte, rte, AP, 5, bwasw_sc_mat);
                free(_cigar);
            } else {
                // for DUP
                s_tlen += 2*(hash_len-dis);
                s_tseq = (uint8_t*)malloc(s_tlen * sizeof(uint8_t));
                ref_offset = at1.offset + APP.seed_len + at1.len_dif + dis - hash_len;
                pac2fa_core(bns, pac, at1.chr, ref_offset-1, &s_tlen, s_tseq);
				s_tlen = s_qlen + dis;
                if (s_tlen < hash_len) {
                    ksw_bi_extend(s_qlen, s_qseq, s_tlen, s_tseq+hash_len-dis, 5, bwasw_sc_mat, abs(s_tlen-s_qlen)+3, hash_len*bwasw_sc_mat[0], hash_len*bwasw_sc_mat[0], AP, &_cigar, &_cigar_n, &_cigar_m);
                    _push_cigar(&s_cigar, &s_clen, &s_cm, _cigar, _cigar_n);
                    free(_cigar);
                }
                else {
                    res = split_indel_map(&s_cigar, &s_clen, &s_cm, s_qseq, s_qlen, s_tseq+hash_len-dis, s_tlen, -dis, AP, hash_num, hash_node);
                }
            }
        } else { //for MIS-MATCH
            s_tlen = s_qlen + dis;
            s_tseq = (uint8_t*)malloc(s_tlen * sizeof(uint8_t));
            ref_offset = at1.offset + APP.seed_len + at1.len_dif;
            pac2fa_core(bns, pac, at1.chr, ref_offset-1, &s_tlen, s_tseq);
			res = ksw_bi_extend(s_qlen, s_qseq, s_tlen, s_tseq, 5, bwasw_sc_mat, abs(s_tlen-s_qlen)+3, 100, 100, AP, &_cigar, &_cigar_n, &_cigar_m);
			_push_cigar(&s_cigar, &s_clen, &s_cm, _cigar, _cigar_n);
			free(_cigar);
        }
    } 
	merge_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), &(la->res[la->cur_res_n].refend), &(la->res[la->cur_res_n].readend), at1.chr, s_cigar, s_clen, s_tlen, s_qlen, bns, pac, read_bseq, AP);
    free(s_qseq); free(s_tseq);
	free(s_cigar);
}

void check_cigar(cigar32_t *cigar, int cigar_len, char *read_name, int read_len)
{
    int i;
    int mid[5] = {0, 0, 0, 0, 0};
    for (i = 0; i < cigar_len; ++i)
        mid[(int)(cigar[i]&0xf)] += (cigar[i]>>4);
    //fprintf(stdout, "\t%d M, %d I, %d D, %d S\t", mid[0], mid[1], mid[2], mid[4]);
    if (mid[0] + mid[1] + mid[4] != read_len) { fprintf(stderr, "\n[cigar len error]: %s cigar error[%d]\t", read_name, read_len); printcigar(stderr, cigar, cigar_len); }
}

int frag_head_bound_fix(frag_msg *f_msg, aln_msg *a_msg, 
                        uint64_t *offset, bntseq_t *bns, uint8_t *pac, 
                        uint8_t *read_bseq, 
                        lsat_aln_per_para APP, lsat_aln_para AP,
                        line_aln_res *la)//aln_res *a_res)
{
    int left_bound = f_msg->frag_left_bound;
    int i, frag_i, seed_x, seed_i, aln_i;
    int read_len, ref_len, read_start/*0-base*/;
    int64_t ref_start;
    if (f_msg->fa_msg[0].srand == 1) {	//'+' srand
        frag_i = f_msg->frag_num-1;
        seed_x = f_msg->fa_msg[frag_i].seed_num - 1;
        seed_i = f_msg->fa_msg[frag_i].seed_i[seed_x];
        aln_i = f_msg->fa_msg[frag_i].seed_aln_i[seed_x];
        if (a_msg[seed_i].read_id != 1) {
            read_len = ((left_bound==0)?0:APP.seed_inv) + (a_msg[seed_i].read_id - left_bound - 1) * APP.seed_step;
            if (read_len < 0) { 
                fprintf(stderr, "[frag_head_bound] Error: %s left bound error.\n", APP.read_name); exit(1); 
            }
            read_start = (left_bound==0)?0:(left_bound * APP.seed_step - APP.seed_inv);
        } else {
            *offset = a_msg[seed_i].at[aln_i].offset;
            la->res[la->cur_res_n].offset = a_msg[seed_i].at[aln_i].offset;
            la->res[la->cur_res_n].refend = la->res[la->cur_res_n].offset - 1;
            la->res[la->cur_res_n].cigar_len = 0;
            return 0;
        }
    } else {	//'-' srand
        frag_i = 0; seed_x = 0;
        seed_i = f_msg->fa_msg[0].seed_i[0];
        aln_i = f_msg->fa_msg[0].seed_aln_i[0];
        read_len = ((left_bound==0)?APP.last_len:APP.seed_inv) + (a_msg[seed_i].read_id - 1 - left_bound) * APP.seed_step;
		if (read_len == 0) {
            *offset = a_msg[seed_i].at[aln_i].offset;
            la->res[la->cur_res_n].offset = a_msg[seed_i].at[aln_i].offset;
            la->res[la->cur_res_n].refend = la->res[la->cur_res_n].offset - 1;
            la->res[la->cur_res_n].cigar_len = 0;
			return 0;
		}
        if (read_len < 0) { 
            fprintf(stderr, "[frag_head_bound] Error: %s left bound error.\n", APP.read_name); exit(1); 
        }
        read_start= (left_bound==0)?0:(APP.last_len+left_bound*APP.seed_step-APP.seed_inv);
    }
    //hash map
    {	
        int hash_len = AP.hash_len, hash_step = AP.hash_step;
        (*offset) = a_msg[seed_i].at[aln_i].offset;
        la->res[la->cur_res_n].offset = a_msg[seed_i].at[aln_i].offset;
        //read
        uint8_t *bseq1 = (uint8_t*)calloc(read_len, sizeof(uint8_t));
        for (i = 0; i < read_len; ++i) bseq1[i] = read_bseq[read_start+i];
        //ref
        ref_len = read_len + hash_step * 2;	//XXX
        ref_start = a_msg[seed_i].at[aln_i].offset - ref_len;	//1-base
        if (ref_start < 1) {
            ref_start = 1;
            ref_len = a_msg[seed_i].at[aln_i].offset - 1;
        }
		uint8_t *bseq2 = (uint8_t*)calloc(ref_len, sizeof(uint8_t));
        pac2fa_core(bns, pac, a_msg[seed_i].at[aln_i].chr, ref_start-1/*0-base*/, &ref_len, bseq2);
		cigar32_t *cigar_=0; int cigar_n_, cigar_m_;
		int qre, tre;
		int res = ksw_extend_r(read_len, bseq1, ref_len , bseq2, 5, bwasw_sc_mat, hash_len, APP.seed_len * bwasw_sc_mat[0], AP, &qre, &tre, &cigar_, &cigar_n_, &cigar_m_);
		if (res != 0) { // not-to-end
			// push head 'S'
			_push_cigar1(&cigar_, &cigar_n_, &cigar_m_, ((read_len-qre)<<4)|CSOFT_CLIP);
		}
		_invert_cigar(&cigar_, cigar_n_);
		(*offset) -= refInCigar(cigar_, cigar_n_);
		(la->res[la->cur_res_n].offset) -= refInCigar(cigar_, cigar_n_);
        la->res[la->cur_res_n].refend = la->res[la->cur_res_n].offset - 1;
        //merge_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), &(la->res[la->cur_res_n].refend), &(la->res[la->cur_res_n].readend), a_msg[seed_i].at[aln_i].chr, *cigar, *cigar_len, refInCigar(*cigar, *cigar_len), readInCigar(*cigar, *cigar_len), bns, pac, read_seq);
		_push_cigar_e(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), &(la->res[la->cur_res_n].refend), &(la->res[la->cur_res_n].readend), cigar_, cigar_n_);
		free(cigar_); free(bseq1); free(bseq2);
    }
    return 0;
}

int frag_tail_bound_fix(frag_msg *f_msg, aln_msg *a_msg, bntseq_t *bns, uint8_t *pac, 
                        uint8_t *read_bseq,
                        lsat_aln_per_para APP, lsat_aln_para AP,
                        line_aln_res *la)
{
    int right_bound = f_msg->frag_right_bound;
    int i, frag_i, seed_x, seed_i, aln_i;
    int read_len, read_start, ref_len;
    int64_t ref_start;

    if (f_msg->fa_msg[0].srand == 1) {	//'+' srand
        seed_i = f_msg->fa_msg[0].seed_i[0];
        aln_i = f_msg->fa_msg[0].seed_aln_i[0];
        read_start = a_msg[seed_i].read_id * APP.seed_step - APP.seed_inv;
        read_len = ((right_bound==APP.seed_all+1)?APP.last_len:APP.seed_inv) + (right_bound - 1 - a_msg[seed_i].read_id) * APP.seed_step;
		if (read_len == 0) return 0;
        if (read_len < 0) { fprintf(stderr, "[frag_tail_bound] Error: %s right bound error.\n", APP.read_name); exit(1); }
    } else {
        frag_i = f_msg->frag_num - 1;
        seed_x = f_msg->fa_msg[frag_i].seed_num - 1;
        seed_i = f_msg->fa_msg[frag_i].seed_i[seed_x];
        aln_i = f_msg->fa_msg[frag_i].seed_aln_i[seed_x];
        if (a_msg[seed_i].read_id != APP.seed_all) {
            read_start = a_msg[seed_i].read_id * APP.seed_step - APP.seed_inv + APP.last_len;
            read_len = ((right_bound== APP.seed_all+1)?0:APP.seed_inv) + (right_bound - 1 - a_msg[seed_i].read_id) * APP.seed_step;
            if (read_len < 0) { fprintf(stderr, "[frag_tail_bound] Error: %s right bound error.\n", APP.read_name); exit(1); }
        } else return 0;
    }
    //hash map
    {
        int hash_len = AP.hash_len, hash_step = AP.hash_step;
        ref_len = read_len + hash_step * 2;
        //read
        uint8_t *bseq1 = (uint8_t*)calloc(read_len, sizeof(uint8_t));
        for (i = 0; i < read_len; ++i)
            bseq1[i] = read_bseq[read_start+i];
        //ref
        ref_start = a_msg[seed_i].at[aln_i].offset + APP.seed_len+ a_msg[seed_i].at[aln_i].len_dif;	//1-base
        uint8_t *bseq2 = (uint8_t*)calloc(ref_len, sizeof(uint8_t));
        pac2fa_core(bns, pac, a_msg[seed_i].at[aln_i].chr, ref_start-1/*0-base*/, &ref_len, bseq2);

		cigar32_t *cigar_=0; int cigar_n_, cigar_m_;
		int qle, tle;
        int res = ksw_extend_c(read_len, bseq1, ref_len, bseq2, 5, bwasw_sc_mat, hash_len, APP.seed_len*bwasw_sc_mat[0], AP, &qle, &tle, &cigar_, &cigar_n_, &cigar_m_);
        if (res != 0) { // not-to-end 
			_push_cigar1(&cigar_, &cigar_n_, &cigar_m_, ((read_len-qle)<<4)|CSOFT_CLIP);
		}
		merge_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), &(la->res[la->cur_res_n].refend), &(la->res[la->cur_res_n].readend), a_msg[seed_i].at[aln_i].chr, cigar_, cigar_n_, refInCigar(cigar_, cigar_n_), readInCigar(cigar_, cigar_n_), bns, pac, read_bseq, AP);
		free(cigar_); free(bseq1); free(bseq2);
    }
    return 0;
}

// generate split alignment, if necessary
// 1. fix split_i : tail_S
// 2. create split_i+1: offset, cigar, head_S, push_cigar
void lsat_res_split(line_aln_res *la, int read_len, lsat_aln_para AP)
{
    int i, j, res_n;
    int op, len, cigar_len; cigar32_t *cigar=0; 
    uint64_t offset; int tail_s, head_s, len1; 

    cigar_len = la->res[0].cigar_len;
    cigar = (cigar32_t*)realloc(cigar, cigar_len * sizeof(cigar32_t));
    for (i = 0; i < cigar_len; ++i) cigar[i] = la->res[0].cigar[i];
    offset = la->res[0].offset;
    res_n = 0, la->res[0].cigar_len = 0;
    for (j = 0; j < cigar_len; ++j) {
        op = cigar[j] & 0xf; 
        if (op == CMATCH) {
            _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), cigar[j]);
        } else if (op == CINS) {
            if ((len = cigar[j] >> 4) >= AP.split_len) { // ins split
                // tail S
                len1 = readInCigar(la->res[res_n].cigar, la->res[res_n].cigar_len);
                tail_s = read_len - len1;
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((tail_s << 4) | CSOFT_CLIP));
                // res_n++
                push_res(la); ++res_n;
                // offset
                la->res[res_n].offset = la->res[res_n-1].offset + refInCigar(la->res[res_n-1].cigar, la->res[res_n-1].cigar_len);
                // head S
                head_s = len + len1;
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((head_s << 4) | CSOFT_CLIP));
            } else _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), cigar[j]);
        } else if (op == CDEL) {
            if ((len = cigar[j] >> 4) >= AP.split_len) { // del split
                // tail S
                len1 = readInCigar(la->res[res_n].cigar, la->res[res_n].cigar_len);
                tail_s = read_len - len1;
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((tail_s << 4) | CSOFT_CLIP));
                // ++res_n
                ++res_n; push_res(la); 
                // offset
                la->res[res_n].offset = la->res[res_n-1].offset + refInCigar(la->res[res_n-1].cigar, la->res[res_n-1].cigar_len) + len;
                // head S
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((len1 << 4) | CSOFT_CLIP));
            } else _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), cigar[j]);
        } else if (op == CSOFT_CLIP) {
            if (j > 0 && j < cigar_len-1 && (cigar[j+1] & 0xf) == CHARD_CLIP) {  // nSmH, mis-match split
                int Sn = cigar[j]>>4, Hn = cigar[j+1]>>4;
                // tail S
                len1 = readInCigar(la->res[res_n].cigar, la->res[res_n].cigar_len);
                tail_s = read_len - len1;
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((tail_s << 4) | CSOFT_CLIP));
                // ++res_n
                ++res_n; push_res(la);
                // offset
                la->res[res_n].offset = la->res[res_n-1].offset + refInCigar(la->res[res_n-1].cigar, la->res[res_n-1].cigar_len) + Hn;
                // head S
                len = Sn;
                head_s = len1 + len;
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((head_s << 4) | CSOFT_CLIP));
				j+=1;
            } else _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), cigar[j]);
        } else if (op != CHARD_CLIP) {
            fprintf(stderr, "\n%s\n[lsat_res_split] Error: Unexpected cigar operation: %d.\n", READ_NAME, op); exit(1);
        }
    }
    free(cigar);
}

void copy_res(res_t *f, res_t *t)
{
    t->offset = f->offset;
    t->chr = f->chr;
    t->nsrand = f->nsrand;
    t->cigar_len = 0;
    _push_cigar(&(t->cigar), &(t->cigar_len), &(t->c_m), f->cigar, f->cigar_len);
    t->refend = f->refend;
    t->readend = f->readend;
    
    t->score = f->score;
    t->mapq = f->mapq;
    t->NM = f->NM;
}

// calculate NM/MAPQ/AS tags
void lsat_res_aux(line_aln_res *la, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, int read_len, lsat_aln_para AP, lsat_aln_per_para APP)
{
    //int read_len = a_res->read_len;
    uint8_t *ref_seq=0; int ref_len;
    int m, i, j, op, len, ref_i, read_i;
    int n_mm, n_m, n_o, n_e; // mis-match/match/gap-open/gap-extension
    int mm_tmp;
    res_t *r; cigar32_t *cigar;

    for (m = 0; m <= la->cur_res_n; ++m) {
        r = la->res+m; 
        cigar = r->cigar;
        ref_len = refInCigar(cigar, r->cigar_len);
        ref_seq = (uint8_t*)realloc(ref_seq, ref_len * sizeof(uint8_t));
        pac2fa_core(bns, pac, r->chr, r->offset-1/*0-base*/, &ref_len, ref_seq);
        ref_i = read_i = n_mm = n_m = n_o = n_e = 0;
        for (i = 0; i < r->cigar_len; ++i) {
            op = cigar[i] & 0xf, len = cigar[i]>>4;
            if (op == CMATCH) {
                mm_tmp = 0;
                for (j = 0; j < len; ++j) {
                    if (read_bseq[read_i++] != ref_seq[ref_i++])
                        ++mm_tmp;
                }
                n_m += (len-mm_tmp);
                n_mm += mm_tmp;
            } else if (op == CINS) {
                read_i+=len;
                n_e += len;
                ++n_o;
            } else if (op == CDEL) {
                ref_i += len;
                n_e += len;
                ++n_o;
            } else if (op == CSOFT_CLIP) {
                read_i+= len;
                continue;
            } else { 
				fprintf(stderr, "\n%s\n[lsat_gen_aux] Error: Unexpected cigar operation: %d.\n", READ_NAME, op); printcigar(stderr, cigar, r->cigar_len); exit(1); 
			}
        }
        if (read_i != read_len || ref_i != ref_len) { 
			fprintf(stderr, "[lsat_gen_aux] Error: %s Unmatched length: read: %d %d.\tref: %d %d\n", APP.read_name, read_i, read_len, ref_i, ref_len); printcigar(stderr, cigar, r->cigar_len); fprintf(stderr, "\n"); exit(1); }
        // calculate NM and AS
        r->NM = n_mm + n_e;
        r->score = n_m * AP.match - n_mm * AP.mis - n_o * AP.gapo - n_e * AP.gape;
        r->mapq = 100; //XXX
		if (r->score < 0) { // XXX
            for (i = m+1; i <= la->cur_res_n; ++i) {
                copy_res(la->res+i, la->res+i-1);
            }
            m--;
            la->cur_res_n--;
		} else {
			la->tol_score += r->score;
			la->tol_NM += r->NM;
		}
    }
	if (la->cur_res_n < 0) la->tol_score = -1;
	else la->tol_score -= (la->cur_res_n * AP.split_pen);
    free(ref_seq);
}

//read_seq: char or uint8_t?
void frag_check(aln_msg *a_msg, frag_msg **f_msg, aln_res *a_res,
        bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, uint8_t **read_rbseq,/*char *read_cseq,*/
        lsat_aln_per_para APP, lsat_aln_para AP,
        int line_n, uint32_t **hash_num, uint64_t ***hash_node) 
{
    strcpy(READ_NAME, APP.read_name);
    int i, j;
    int max_len = APP.read_len;
    uint8_t *bseq1 = (uint8_t*)malloc((max_len+1)*sizeof(uint8_t)); // for qurey bseq and target bseq
    uint8_t *bseq2 = (uint8_t*)malloc((max_len+1)*sizeof(uint8_t)); int ref_max_blen = max_len+1;

    // realloc mem, if necessary
    if (line_n > a_res->l_m) {
        extern void aln_reloc_res(aln_res *a_res, int line_n, int XA_m);
        aln_reloc_res(a_res, line_n, AP.res_mul_max);
    }
         
    // initialization
    a_res->l_n = line_n;
    for (i = 0; i < line_n; ++i) {
        a_res->la[i].line_score = (*f_msg)[i].line_score;
        a_res->la[i].cur_res_n = 0;
        a_res->la[i].tol_score = 0; a_res->la[i].tol_NM = 0;
        for (j = 0; j < a_res->la[i].res_m; ++j) a_res->la[i].res[j].cigar_len = 0;
    }

    int res_n;
    uint64_t offset;
    // check for every line
    for (j = 0; j < line_n; ++j) {
        line_aln_res *la = a_res->la+j;
        la->res[la->cur_res_n].cigar_len = 0;
        la->res[la->cur_res_n].nsrand = (((*f_msg)+j)->fa_msg[0].srand == 1)?1:0;
        la->res[la->cur_res_n].chr = ((*f_msg)+j)->fa_msg[0].chr;
        la->res[la->cur_res_n].readend = 0;
        /* extend once
           for every frag : take it as a RIGHT frag
           for every interval : SV breakpoint(s).*/
        if (((*f_msg)+j)->fa_msg[0].srand == 1) {
            //head 'S'
            if (((*f_msg)+j)->frag_left_bound > 0) {
                cigar32_t s_cigar = ((((*f_msg)+j)->frag_left_bound * APP.seed_step - APP.seed_inv) << 4) | CSOFT_CLIP;
                res_n = la->cur_res_n;
                _push_cigar_e1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), 0, &(la->res[res_n].readend), s_cigar);
            }
            //fix the boundary blank before first frag, if it exists
            frag_head_bound_fix((*f_msg)+j, a_msg, &offset, bns, pac, read_bseq, APP, AP, la);
            for (i = ((*f_msg)+j)->frag_num-1; i > 0; --i) {
                frag_extend((*f_msg)+j, a_msg, i, bns, pac, read_bseq, bseq1, &bseq2, &ref_max_blen, APP, AP, la);
                split_mapping(bns, pac, read_bseq, ((*f_msg)+j), a_msg, hash_num, hash_node, i, i-1, APP, AP, la);
            }
            frag_extend((*f_msg)+j, a_msg, i, bns, pac, read_bseq, bseq1, &bseq2, &ref_max_blen, APP, AP, la);
            //fix the boundary blank after the last frag
            frag_tail_bound_fix((*f_msg)+j, a_msg, bns, pac, read_bseq, APP, AP, la);
            //tail 'S'
            if (((*f_msg)+j)->frag_right_bound <= APP.seed_all) {
                cigar32_t s_cigar = ((APP.read_len - (((*f_msg)+j)->frag_right_bound-1) * APP.seed_step) << 4) | CSOFT_CLIP; 
                res_n = la->cur_res_n;
                _push_cigar_e1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), 0, &(la->res[res_n].readend), s_cigar);
            }
            //split and aux
            lsat_res_split(a_res->la+j, APP.read_len, AP);
            lsat_res_aux(a_res->la+j, bns, pac, read_bseq, APP.read_len, AP, APP);
        } else {	//'-' srand
            //convert into rev-com
            if (*read_rbseq == NULL) {
                *read_rbseq = (uint8_t*)calloc(APP.read_len, sizeof(uint8_t));
                for (i = 0; i < APP.read_len; ++i) (*read_rbseq)[i] =  (read_bseq[APP.read_len-1-i]<4)?3-read_bseq[APP.read_len-1-i]:4;
            }
            for (i = 0; i < APP.seed_all; ++i) a_msg[i].read_id = (APP.seed_all + 1 - a_msg[i].read_id);
            int tmp = ((*f_msg)+j)->frag_left_bound;
            ((*f_msg)+j)->frag_left_bound = APP.seed_all + 1 - ((*f_msg)+j)->frag_right_bound;
            ((*f_msg)+j)->frag_right_bound = APP.seed_all + 1 - tmp;
            //head 'S'
            if (((*f_msg)+j)->frag_left_bound > 0) {
                cigar32_t s_cigar = ((((*f_msg)+j)->frag_left_bound * APP.seed_step - APP.seed_inv + APP.last_len) << 4) | CSOFT_CLIP;
                res_n = la->cur_res_n;
                _push_cigar_e1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), 0, &(la->res[res_n].readend), s_cigar);
            }
            //fix the boundary blank before first frag
            frag_head_bound_fix((*f_msg)+j, a_msg, &offset, bns, pac, *read_rbseq, APP, AP, la);
            for (i = 0; i < ((*f_msg)+j)->frag_num-1; ++i) {
                frag_extend((*f_msg)+j, a_msg, i, bns, pac, *read_rbseq, bseq1, &bseq2, &ref_max_blen, APP, AP, la);
                split_mapping(bns, pac, *read_rbseq, (*f_msg)+j, a_msg, hash_num, hash_node, i, i+1, APP, AP, la);
            }
            frag_extend(((*f_msg)+j), a_msg, i, bns, pac, *read_rbseq, bseq1, &bseq2, &ref_max_blen,  APP, AP, la);
            //fix the boundary blank after the last frag
            frag_tail_bound_fix(((*f_msg)+j), a_msg, bns, pac, *read_rbseq, APP, AP, la);
            //tail 'S'
            if (((*f_msg)+j)->frag_right_bound <= APP.seed_all) {
                cigar32_t s_cigar = (((APP.seed_all - ((*f_msg)+j)->frag_right_bound+1) * APP.seed_step - APP.seed_inv) << 4) | CSOFT_CLIP;
                res_n = la->cur_res_n;
                _push_cigar_e1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), &(la->res[res_n].refend), &(la->res[res_n].readend), s_cigar);
            }
            lsat_res_split(a_res->la+j, APP.read_len, AP);
            lsat_res_aux(a_res->la+j, bns, pac, *read_rbseq, APP.read_len, AP, APP);
            for (i = 0; i < APP.seed_all; ++i) a_msg[i].read_id = (APP.seed_all + 1 - a_msg[i].read_id);
        }
    }
    // filter line-aln-res in merged-lines, get opt-res OR multi-sub-opt-res
    //res_filter(a_res);

    free(bseq1); free(bseq2);
}
