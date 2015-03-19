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
#include "ssw.h"
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
	for (i = 0; i < tf_msg->frag_max; ++i)
	{
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
		for (j = 0; j < f_msg[i].frag_max; ++j)
		{
			free(f_msg[i].fa_msg[j].cigar);
			free(f_msg[i].fa_msg[j].seed_i);
			free(f_msg[i].fa_msg[j].seed_aln_i);
		}
		free(f_msg[i].fa_msg);
	}
	free(f_msg);
}

int frag_set_msg(aln_msg *a_msg, int seed_i, int aln_i,
				 int FLAG, frag_msg *f_msg, int frag_i)//FLAG 0:start / 1:end / 2:seed
{
	//if (frag_i == 0)
	//	f_msg->seed_num=0;
	if (FLAG == FRAG_END) {	//end
		f_msg->fa_msg[frag_i].chr = a_msg[seed_i].at[aln_i].chr;
		f_msg->fa_msg[frag_i].srand = a_msg[seed_i].at[aln_i].nsrand;	//+:1/-:-1
		f_msg->fa_msg[frag_i].seed_i[0] = seed_i;
		f_msg->fa_msg[frag_i].seed_aln_i[0] = aln_i;
		f_msg->fa_msg[frag_i].seed_num = 1;
		//f_msg->seed_num++;
		f_msg->fa_msg[frag_i].flag = UNCOVERED;
		f_msg->fa_msg[frag_i].cigar_len = 0;
        //trg_n
        f_msg->fa_msg[frag_i].trg_n = 0;
	}
	else if(FLAG==FRAG_START) {	//start
		f_msg->frag_num = frag_i + 1;
	}
	else {			//seed
		f_msg->fa_msg[frag_i].seed_i[f_msg->fa_msg[frag_i].seed_num] = seed_i;
		f_msg->fa_msg[frag_i].seed_aln_i[f_msg->fa_msg[frag_i].seed_num] = aln_i;
		//f_msg->seed_num++;
		f_msg->fa_msg[frag_i].seed_num++;
	}
	return 0;
}

int frag_trg_set(frag_dp_node f_node, frag_msg *f_msg, int frag_i)
{
    f_msg->fa_msg[frag_i].trg_n |= f_node.trg_n;
    if (f_node.trg_n & 0x1)
        f_msg->fa_msg[frag_i].next_trg = f_node.next_trg;
    else // 0x2
        f_msg->fa_msg[frag_i].pre_trg = f_node.pre_trg;
    return 0;
}

//XXX combine get_ref_intv & get_read_intv
int get_ref_intv(uint8_t *seq1, int *N_flag, 
				 bntseq_t *bns, uint8_t *pac, 
				 aln_msg *a_msg, 
				 int seed1_i, int seed1_aln_i, int seed2_i, int seed2_aln_i, 
                 lsat_aln_per_para *APP)
{
	int64_t start; int32_t len; int N_len;
    start = a_msg[seed1_i].at[seed1_aln_i].offset + APP->seed_len - 1 + a_msg[seed1_i].at[seed1_aln_i].len_dif;
    len = a_msg[seed2_i].at[seed2_aln_i].offset - 1 - start;
    if (len <= 0) return 0;
    //printf("offset: %lld ", (long long)start);
    pac2fa_core(bns, pac, a_msg[seed1_i].at[seed1_aln_i].chr, start, &len, 1, N_flag, &N_len, seq1);
    //int j;
	//fprintf(stdout, "ref:\t%d\t%ld\t", len, start);
	//for (j = 0; j < len; j++)
	//	fprintf(stdout, "%d", (int)seq1[j]);
	//fprintf(stdout, "\n");
	return (int)len;
}

//@func: convert '0123' to 'ACGT', check for 'N'
int get_read_seq(uint8_t *read_seq, char *read_char, int start/*0-base*/, int *len, int *flag)
{
	int i;
	for (i = 0; i < *len; ++i)
	{
		read_char[i] = nst_nt4_table[(int)read_seq[start+i]];
	}
    return 0;
}

int get_read_intv(uint8_t *seq2, char *read_seq, 
				  aln_msg *a_msg, 
				  int seed1_i, int seed1_aln_i, int seed2_i, int seed2_aln_i, 
				  int *band_width, 
				  lsat_aln_per_para *APP)
{
	int32_t i;
	int j;

	//set band-width
	*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
	
	
	//convert char seq to int seq

	if (a_msg[seed1_i].at[seed1_aln_i].nsrand == 1)
	{
		*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
		for (j=0, i = a_msg[seed1_i].read_id * APP->seed_step - APP->seed_inv; i < (a_msg[seed2_i].read_id-1) * APP->seed_step; ++j, ++i)
			seq2[j] = nst_nt4_table[(int)read_seq[i]];
	}
	else
	{
		//printf("read: %d, %d ", a_msg[seed1_i].read_id,a_msg[seed1_i].read_id * 2 * seed_len);
		*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
        for (j=0, i = APP->last_len+a_msg[seed1_i].read_id*APP->seed_step-APP->seed_inv; i < APP->last_len+(a_msg[seed2_i].read_id-1)*APP->seed_step; ++j, ++i)
		//for (j=0, i = a_msg[seed1_i].read_id * APP->seed_step; i < (a_msg[seed2_i].read_id * 2 - 1) * seed_len; ++j, ++i)
		{
			seq2[j] = nst_nt4_table[(int)read_seq[i]];
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

    //fprintf(stdout, "cigar %d: ", cigar_len);
	for (i = 0; i < cigar_len; i++)
		fprintf(outp, "%d%c", (int)(cigar[i]>>4), "MIDNSHP=X"[(int)(cigar[i] & 0xf)]);
	//fprintf(stdout, "\t");
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

//merge adjacent same operations XXX
void _push_cigar(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t *_cigar, int _cigar_n)
{
    if (_cigar_n == 0) return;
	int i,j;
	i = *cigar_n;
	if (((i-1) >= 0) && (((*cigar)[i-1] & 0xf) == (_cigar[0] & 0xf)))
	{
		(*cigar)[i-1] += ((_cigar[0] >> 4) << 4);
		j = 1;
	}
	else j = 0;

	for (; j < _cigar_n; ++i,++j) {
		if (i == *cigar_m) {
			(*cigar_m) <<= 1;
			(*cigar) = (cigar32_t*)realloc(*cigar, (*cigar_m) * sizeof (cigar32_t));
			if ((*cigar) == NULL)	{fprintf(stderr, "\n[frag_check] Memory is not enougy.\n");exit(-1);}
		}
		(*cigar)[i] = _cigar[j];
	}
	*cigar_n = i;
}

void _push_cigar1(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t _cigar) {
	if (_cigar >> 4 == 0) return;
    int i;
    i = *cigar_n;
    if (((i-1) >=0) && (((*cigar)[i-1] & 0xf) == (_cigar & 0xf)))
        (*cigar)[i-1] += ((_cigar >> 4) << 4);
    else
    {
        if (i == *cigar_m) {
            (*cigar_m) <<= 1;
			(*cigar) = (cigar32_t*)realloc(*cigar, (*cigar_m) * sizeof (cigar32_t));
			if ((*cigar) == NULL)	{fprintf(stderr, "\n[frag_check] Memory is not enougy.\n");exit(-1);}
        }
        (*cigar)[i] = _cigar;
        ++(*cigar_n);
    }
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
            if (i < 0) { fprintf(stderr, "\n[pop_cigar] Cigar error,\n"); exit(0); }
        }
    }
    else 
    {
        fprintf(stderr, "\n[pop_cigar] Unknown pop mode. (%d)\n", pop_mode);
        exit(0);
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

//return 'MIS' length in cigar
int readInCigar(cigar32_t *cigar, int cigar_len) {
    int i, read_len=0;
    for (i = 0; i < cigar_len; ++i) {
        if ((cigar[i] & 0xf) == CMATCH || (cigar[i] & 0xf) == CINS || (cigar[i] & 0xf) == CSOFT_CLIP)
            read_len += (cigar[i]  >> 4);
        else if ((cigar[i] & 0xf) != CDEL && (cigar[i] & 0xf) != CHARD_CLIP) {
            fprintf(stderr, "\n[readInCigar] Cigar Error.\n");
            printcigar(stderr, cigar, cigar_len); 
            exit(0);
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
            fprintf(stderr, "\n[readInCigar] Cigar Error.\n");
            printcigar(stderr, cigar, cigar_len); 
            exit(0);
        }
    }
    return read_len;
}

int cut_cigar(cigar32_t **cigar, int *cigar_len, line_aln_res *la, int read_len) {
    // cut aln-cigar into split-cigar on 'nSmH' 
    // TODO: multi-'nSmH', S len, '+/-' strand 
    int i;
    cigar32_t pre_soft_cigar, next_soft_cigar;
    int pre_clen, next_clen, cut_flag;
    int next_soft_len, pre_soft_len;
    uint64_t next_offset, tmp_offset;

    // printcigar(stdout, *cigar, *cigar_len);
    while (1) { // for multi-'nSmH'
        pre_clen = 0, next_clen = 0, cut_flag = 0;
        next_soft_len = 0, pre_soft_len = 0;
        next_soft_len += readInCigar(la->res[la->cur_res_n].cigar, la->res[la->cur_res_n].cigar_len);
        tmp_offset = la->res[la->cur_res_n].offset+readInCigar(la->res[la->cur_res_n].cigar, la->res[la->cur_res_n].cigar_len);
        for (i = 0; i < *cigar_len; ++i) {
            if ((i < *cigar_len-1) && (((*cigar)[i]) & 0xf) == CSOFT_CLIP && (((*cigar)[i+1]) & 0xf) == CHARD_CLIP) {
                // pre
                pre_clen = i+1;
                _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), 
                        &(la->res[la->cur_res_n].c_m), *cigar, pre_clen);
                // next
                next_offset = tmp_offset;
                next_offset += (((*cigar)[i+1]) >> 4);
                next_soft_len += (((*cigar)[i]) >> 4);
                cut_flag = 1;
                _deQueue_cigar(cigar, cigar_len, pre_clen+1);
                break;
            } else if ((((*cigar)[i]) & 0xf) == CMATCH) {
                tmp_offset += (((*cigar)[i]) >> 4);
                next_soft_len += (((*cigar)[i]) >> 4);
            } else if ((((*cigar)[i]) & 0xf) == CSOFT_CLIP || (((*cigar)[i]) & 0xf) == CINS) {
                next_soft_len += (((*cigar)[i]) >> 4);
            } else if ((((*cigar)[i]) & 0xf) ==  CDEL) {
                tmp_offset += (((*cigar)[i]) >> 4);
            } else {
                fprintf(stderr, "\n[split_mapping] Cigar Error.\n");
                printcigar(stderr, *cigar, *cigar_len);
                exit(0);
            }
        }
        if (cut_flag) {
            // pre
            pre_soft_len = read_len - readInCigar(la->res[la->cur_res_n].cigar, la->res[la->cur_res_n].cigar_len);
            pre_soft_cigar = ((pre_soft_len << 4) | CSOFT_CLIP);
            _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), &pre_soft_cigar, 1);
            // next
            ++(la->cur_res_n);
            if (la->cur_res_n == la->res_m) {
                la->res_m <<= 1;
                la->res = (res_t*)realloc(la->res, la->res_m * sizeof(res_t));
                for (i = la->res_m >> 1; i < la->res_m; ++i) {
                    la->res[i].c_m = 100;
                    la->res[i].cigar = (cigar32_t*)malloc(100 * sizeof(cigar32_t));
                }
            }
            la->res[la->cur_res_n].offset = next_offset;
            la->res[la->cur_res_n].chr = la->res[la->cur_res_n-1].chr;
            la->res[la->cur_res_n].nsrand = la->res[la->cur_res_n-1].nsrand;
            la->res[la->cur_res_n].cigar_len = 0;
            next_soft_cigar = (next_soft_len << 4 | CSOFT_CLIP);
            _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), &next_soft_cigar, 1);
        } else break;
    }
    return 0;
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
//Before, frag has the first seed's msg already
void merge_cigar(frag_msg *f_msg, int f_i, 
				 cigar32_t *aln_cigar, int cigar_len, 
				 int reflen, int readlen, 
				 bntseq_t *bns, uint8_t *pac, 
				 char *read_seq, int read_len)
{
    if (cigar_len == 0) return;
	cigar32_t *cigar = (cigar32_t*)malloc(cigar_len * sizeof(cigar32_t));
	int ii;
	for (ii= 0; ii < cigar_len; ++ii) cigar[ii] = aln_cigar[ii];

	int *fc_len = &(f_msg->fa_msg[f_i].cigar_len);
	cigar32_t **fcigar = &(f_msg->fa_msg[f_i].cigar);
	//refresh cigar msg
	//XXX
	if (((*fcigar)[*fc_len-1] & 0xf) != CMATCH || (cigar[0] & 0xf) != CMATCH)	//bound-repair
	{
		    /* seq1, ref */ /*    seq2, read     */
		int len1, len11=0,  len2, len21=0, len22=0, len_dif1=0, len_dif2=0;
		int bd_cigar_len=0, b=0, min_b, score, ci=0;
		int N_flag, N_len;
		int i, j, left=1, right=1;
		uint8_t *seq1, *seq2;
		cigar32_t *bd_cigar = 0;
		int64_t ref_start;
		int32_t read_start;
		//XXX 
		int md = 3;

		while (1) {	//get a right boundary-cigar(without 'I'/'D' on the head or tail)
			//extend 3-match boundary
			if (left) {
				while ((*fc_len) >= 1) {	// get a nM, where n > 3
					if (((*fcigar)[*fc_len-1] & 0xf) == CMATCH && ((*fcigar)[*fc_len-1]>>4) > md) {
						(*fcigar)[*fc_len-1] -= (md << 4); len21 += md; break;}
					else if (((*fcigar)[*fc_len-1] & 0xf) == CMATCH) {
						len21 += (int)((*fcigar)[*fc_len-1] >> 4); --(*fc_len); }
					else if (((*fcigar)[*fc_len-1] & 0xf) == CINS) {
						len21 += (int)((*fcigar)[*fc_len-1] >> 4); len_dif1 -= (int)((*fcigar)[*fc_len-1] >> 4); b += (int)((*fcigar)[*fc_len-1]>>4); --(*fc_len); }
					else if (((*fcigar)[*fc_len-1] & 0xf) == CDEL) {
						len_dif1 += (int)((*fcigar)[*fc_len-1] >> 4); b += (int)((*fcigar)[*fc_len-1] >> 4); --(*fc_len); }
					else { fprintf(stderr, "\n[merge cigar] Unexpected(1): "); printcigar(stderr, (*fcigar), (*fc_len)); exit(-1);}
				}
				len11 = len21 + len_dif1;
				read_start = f_msg->fa_msg[f_i].cigar_read_end - len21 + 1;	//1-based
				ref_start = f_msg->fa_msg[f_i].cigar_ref_end - len11 + 1;	//1-based
				//printf("read_start: %d ref_start: %lld ", read_start, (long long)ref_start);
			}
			//extend 3-match boundary
			if (right) {
				while (ci < cigar_len) {	// get a nM, where n > 3
					if ((cigar[ci] & 0xf) == CMATCH && (cigar[ci] >> 4) > md) {
						cigar[ci] -= (md << 4); len22 += md; break;}
					else if ((cigar[ci] & 0xf) == CMATCH) {
						len22 += (int)(cigar[ci] >> 4); ci++; }
					else if ((cigar[ci] & 0xf) == CINS) {
						len22 += (int)(cigar[ci] >> 4); len_dif2 -= (int)(cigar[ci]>>4); b += (int)(cigar[ci]>>4); ++ci; }
					else if ((cigar[ci] & 0xf) == CDEL) {
						len_dif2 += (int)(cigar[ci] >> 4); b += (int)(cigar[ci]>>4); ++ci; }
					else { fprintf(stderr, "\n[merge cigar] Unexpected(2): "); printcigar(stderr, (cigar+ci), cigar_len-ci); exit(-1);}
				}
			}
			len2 = len21+len22; len1 = len2 + len_dif1+len_dif2;
			//printf("len1: %d len21: %d len22: %d len_dif %d len2: %d\t", len1, len21, len22, len_dif, len2);
			min_b = abs(len_dif1+len_dif2) + md; b = MAXOFTWO(b, min_b);
			seq1 = (uint8_t *)malloc((len1+1) * sizeof(uint8_t)); seq2 = (uint8_t *)malloc((len2+1) * sizeof(uint8_t));
			pac2fa_core(bns, pac, f_msg->fa_msg[f_i].chr, ref_start-1, &len1, 1, &N_flag, &N_len, seq1);
			for (i=0, j=read_start-1; j<read_start+len2-1; ++i, ++j) seq2[i] = nst_nt4_table[(int)read_seq[j]];
			//XXX score check
			score = ksw_global(len2, seq2, len1, seq1, 5, bwasw_sc_mat, 5, 2, b, &bd_cigar_len, &bd_cigar);
			//printf("c");printcigar(bd_cigar, bd_cigar_len);
			free(seq1); free(seq2);
			if ((bd_cigar[0] & 0xf) == CMATCH) left = 0;
			if ((bd_cigar[bd_cigar_len-1] & 0xf) == CMATCH) right = 0;

			if ((left + right) == 0 )
			{
				//printf("left: %d right %d\n", left, right);
				break;
			}
			if (ci == cigar_len) 
			{
				//printf("ci: %d\n", ci); 
				break;
			} 
			left = right = 1;
			free(bd_cigar);
		}
		_push_cigar(fcigar, fc_len, &(f_msg->fa_msg[f_i].cigar_max), bd_cigar, bd_cigar_len);
		_push_cigar(fcigar, fc_len, &(f_msg->fa_msg[f_i].cigar_max), cigar+ci, cigar_len-ci);
		free(bd_cigar);
	}
	else _push_cigar(fcigar, fc_len, &(f_msg->fa_msg[f_i].cigar_max), cigar, cigar_len);
	//refresh coordinates msg
	f_msg->fa_msg[f_i].cigar_ref_end += reflen;
	f_msg->fa_msg[f_i].cigar_read_end += readlen;
	f_msg->fa_msg[f_i].len_dif += (reflen-readlen);
	free(cigar);
}

void push_trg(line_node **trg, int *trg_n, int *trg_m, line_node t)
{
    if (*trg_n == *trg_m) {
        (*trg_m) <<= 1;
        (*trg) = (line_node*)realloc(*trg, (*trg_m) * sizeof(line_node));
        if (*trg == NULL) { fprintf(stderr, "[push trg] Error: Not enough memory.\n"); exit(0); }
    }
    (*trg)[(*trg_n)++] = t;
}
//Extend the intervals between seeds in the same frag
//XXX un-seed region: end of + srand and start of - end 
int frag_extend(frag_msg *f_msg, aln_msg *a_msg, int f_i, 
				bntseq_t *bns, uint8_t *pac, char *read_seq, int read_len,
				uint8_t *seq1, uint8_t *seq2, 
                lsat_aln_per_para *APP, lsat_aln_para *AP,
				line_aln_res *la)
{
	int i, seed_i, seed_aln_i, last_i, last_aln_i;
	int len1, len2;
	int N_flag;
	int b, min_b, cigar_len, score;

	//banded global alignment for intervals between seeds
	cigar_len = 0;
	
	//get the first node
	if (f_msg->fa_msg[f_i].srand == 1)
	{
		i = f_msg->fa_msg[f_i].seed_num-1;
		last_i = f_msg->fa_msg[f_i].seed_i[i];
		last_aln_i = f_msg->fa_msg[f_i].seed_aln_i[i];
        f_msg->fa_msg[f_i].cigar_read_start = (a_msg[last_i].read_id-1) * APP->seed_step + 1;
        f_msg->fa_msg[f_i].cigar_read_end = f_msg->fa_msg[f_i].cigar_read_start - 1 + APP->seed_len;
	} else {	//'-' srand 
		last_i = f_msg->fa_msg[f_i].seed_i[0];
		last_aln_i = f_msg->fa_msg[f_i].seed_aln_i[0];
        f_msg->fa_msg[f_i].cigar_read_start = APP->last_len + (a_msg[last_i].read_id-1) * APP->seed_step + 1; 
        f_msg->fa_msg[f_i].cigar_read_end = f_msg->fa_msg[f_i].cigar_read_start - 1 + APP->seed_len;
	}

	//copy the first seed's cigar to frag
	{
		_push_cigar(&(f_msg->fa_msg[f_i].cigar), &(f_msg->fa_msg[f_i].cigar_len), &(f_msg->fa_msg[f_i].cigar_max), a_msg[last_i].at[last_aln_i].cigar, a_msg[last_i].at[last_aln_i].cigar_len);
		f_msg->fa_msg[f_i].len_dif = a_msg[last_i].at[last_aln_i].len_dif;
		//start/end: 1-base
		f_msg->fa_msg[f_i].cigar_ref_start = a_msg[last_i].at[last_aln_i].offset;
		f_msg->fa_msg[f_i].cigar_ref_end = a_msg[last_i].at[last_aln_i].offset + APP->seed_len -1 + a_msg[last_i].at[last_aln_i].len_dif;
	}
	if (f_msg->fa_msg[f_i].srand == 1)
	{
		for (--i; i >= 0; --i)
		{
			seed_i = f_msg->fa_msg[f_i].seed_i[i];
			seed_aln_i = f_msg->fa_msg[f_i].seed_aln_i[i];
			//get ref seq and seed-interval seq
			if ((len1 = get_ref_intv(seq1, &N_flag, bns, pac, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, APP)) < 0) { fprintf(stderr, "\n[frag extend] frag path error : (%d,%d) -> %d %lld %d %lld\n", last_i, last_aln_i, a_msg[last_i].read_id, (long long)a_msg[last_i].at[last_aln_i].offset, a_msg[seed_i].read_id, (long long)a_msg[seed_i].at[seed_aln_i].offset); return -1;}
			len2 = get_read_intv(seq2, read_seq, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, &b, APP);
			min_b = abs(len2-len1)+3;
			b = MAXOFTWO(b, min_b);
			cigar_len = 0;
			cigar32_t *cigar=0;
			//XXX score check
			score = ksw_global(len2, seq2, len1, seq1, 5, bwasw_sc_mat, 5, 2, b, &cigar_len, &cigar);
            //printf("score1: %d\n", score);
			merge_cigar(f_msg, f_i, cigar, cigar_len, len1, len2, bns, pac, read_seq, read_len);
			merge_cigar(f_msg, f_i, a_msg[seed_i].at[seed_aln_i].cigar, a_msg[seed_i].at[seed_aln_i].cigar_len, APP->seed_len+a_msg[seed_i].at[seed_aln_i].len_dif, APP->seed_len, bns, pac, read_seq, read_len);//merge seed to frag
			last_i = seed_i;
			last_aln_i = seed_aln_i;
			free(cigar);
		}
	}
	else 
	{
		for (i = 1; i < f_msg->fa_msg[f_i].seed_num; ++i)
		{
			seed_i = f_msg->fa_msg[f_i].seed_i[i];
			seed_aln_i = f_msg->fa_msg[f_i].seed_aln_i[i];
			//get ref seq and seed-interval seq
			if ((len1 = get_ref_intv(seq1, &N_flag, bns, pac, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, APP)) < 0) { fprintf(stderr, "\n[frag extend] frag path error : (%d,%d) -> %d %lld %d %lld\n", last_i, last_aln_i, a_msg[last_i].read_id, (long long)a_msg[last_i].at[last_aln_i].offset, a_msg[seed_i].read_id, (long long)a_msg[seed_i].at[seed_aln_i].offset); return -1;}
			len2 = get_read_intv(seq2, read_seq, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, &b, APP);
			min_b = abs(len2-len1)+3;
			b = MAXOFTWO(b, min_b);
			cigar_len = 0;
			cigar32_t *cigar=0;
			//XXX score check
			score = ksw_global(len2, seq2, len1, seq1, 5, bwasw_sc_mat, 5, 2, b, &cigar_len, &cigar);
            //printf("score2: %d\n", score);
			//merge_cigar
			merge_cigar(f_msg, f_i, cigar, cigar_len, len1, len2, bns, pac, read_seq, read_len);
			merge_cigar(f_msg, f_i, a_msg[seed_i].at[seed_aln_i].cigar, a_msg[seed_i].at[seed_aln_i].cigar_len, APP->seed_len+a_msg[seed_i].at[seed_aln_i].len_dif, APP->seed_len, bns, pac, read_seq, read_len);//merge seed to frag

			last_i = seed_i;
			last_aln_i = seed_aln_i;
			free(cigar);
		}
	}
    _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), f_msg->fa_msg[f_i].cigar, f_msg->fa_msg[f_i].cigar_len);
	return 0;
}

//for '-' srand:
//1. convert read seq into reverse-complement seq.
//2. split-map the re-co seq to the ref.
//3. combine with aln results of other regions.
//@func: deal with seqs between 2 frags
void split_mapping(cigar32_t **split_cigar, int *split_len, int *split_m,
				   bntseq_t *bns, uint8_t *pac, 
				   char *read_seq, 
				   frag_msg *f_msg, aln_msg *a_msg, 
				   uint32_t **hash_num, uint64_t ***hash_node, 
				   int f1_i, int f2_i, 
                   lsat_aln_per_para *APP, lsat_aln_para *AP,
				   line_aln_res *la)
{
    int i;
    int nsrand, s1_i, s1_aln_i, seed_num, s2_i, s2_aln_i;

	nsrand = f_msg->fa_msg[f1_i].srand;
	if (nsrand == 1) {
		s1_i = f_msg->fa_msg[f1_i].seed_i[0];
		s1_aln_i = f_msg->fa_msg[f1_i].seed_aln_i[0];

		seed_num = f_msg->fa_msg[f2_i].seed_num;
		s2_i = f_msg->fa_msg[f2_i].seed_i[seed_num-1];
		s2_aln_i = f_msg->fa_msg[f2_i].seed_aln_i[seed_num-1];
	} else {
		seed_num = f_msg->fa_msg[f1_i].seed_num;
		s1_i = f_msg->fa_msg[f1_i].seed_i[seed_num-1];
		s1_aln_i = f_msg->fa_msg[f1_i].seed_aln_i[seed_num-1];

		s2_i = f_msg->fa_msg[f2_i].seed_i[0];
		s2_aln_i = f_msg->fa_msg[f2_i].seed_aln_i[0];
	}

	int split_read_len, split_ref_len, res;
	uint8_t *split_read_seq, *split_ref_seq;
	int64_t ref_offset;
	//init value for hash-map
	int hash_len = AP->hash_len, key_len = AP->hash_key_len, hash_step = AP->hash_step;
	int nt_n = NT_N;
	int hash_size = (int)pow(nt_n, key_len);
	int N_flag, N_len;

	split_read_len = (a_msg[s2_i].read_id-a_msg[s1_i].read_id) * APP->seed_step - APP->seed_len;
	split_read_seq = (uint8_t*)malloc(split_read_len * sizeof(uint8_t));

	//if (nsrand == 1) for (i = 0; i < split_read_len; ++i) split_read_seq[i] = nst_nt4_table[(int)read_seq[a_msg[s1_i].read_id*APP->seed_step-APP->seed_inv+i]];
    //else for (i = 0; i < split_read_len; ++i) split_read_seq[i] = nst_nt4_table[(int)read_seq[a_msg[s1_i].read_id*APP->seed_step + i]];//'-' srand	

    int bd;
    get_read_intv(split_read_seq, read_seq, a_msg, s1_i, s1_aln_i, s2_i, s2_aln_i, &bd, APP);

	//check SV-type
	if (a_msg[s1_i].at[s1_aln_i].chr == a_msg[s2_i].at[s2_aln_i].chr &&
		a_msg[s1_i].at[s1_aln_i].nsrand == a_msg[s2_i].at[s2_aln_i].nsrand)
	{
		int64_t pos = a_msg[s1_i].at[s1_aln_i].offset+APP->seed_len-1+a_msg[s1_i].at[s1_aln_i].len_dif;
		int64_t exp = a_msg[s1_i].at[s1_aln_i].offset + (a_msg[s2_i].read_id - a_msg[s1_i].read_id) * APP->seed_step;	
		int64_t act = a_msg[s2_i].at[s2_aln_i].offset;
		int dis = (act-exp) - a_msg[s1_i].at[s1_aln_i].len_dif;

		//fprintf(stdout, "s1: %lld, exp_s2: %lld act_s2: %lld %d\n", (long long)pos, (long long)exp, (long long)act, dis);
		if (dis > APP->match_dis)
		{	
#ifdef __DEBUG__
			fprintf(stdout, "%d\t%lld\t%d\tDEL\n", a_msg[s1_i].at[s1_aln_i].chr, (long long)pos, dis);
#endif
			//long-del XXX
			//if (dis < split_read_len)
			{
				split_ref_len = split_read_len + dis;
				split_ref_seq = (uint8_t*)malloc(split_ref_len * sizeof(uint8_t));
				ref_offset = a_msg[s1_i].at[s1_aln_i].offset + APP->seed_len + a_msg[s1_i].at[s1_aln_i].len_dif;
				//XXX s1/s2???
				pac2fa_core(bns, pac, a_msg[s1_i].at[s1_aln_i].chr, ref_offset-1, &split_ref_len, 1/*+-*/, &N_flag, &N_len, split_ref_seq);

                if (split_read_len < hash_len)
                {
                    cigar32_t *k_cigar=0; int k_clen;
                    int score = ksw_global(split_read_len, split_read_seq, split_ref_len, split_ref_seq, 5, bwasw_sc_mat, 5, 2, dis+3, &k_clen, &k_cigar);
                    *split_len=0;
                    _push_cigar(split_cigar, split_len, split_m, k_cigar, k_clen);
                    free(k_cigar);
                }
                res = split_indel_map(split_cigar, split_len, split_m, split_read_seq, split_read_len, split_ref_seq, split_ref_len, ref_offset, hash_len, hash_step, AP->split_len,  hash_num, hash_node, key_len, hash_size);
                if (res & 1) // pull trigger, cut cigar
                {
                    //printf("del-score: 1\t"); printcigar(stdout, *split_cigar, *split_len); printf("\n");
                    if (nsrand == 1) {
                        if (f_msg->fa_msg[f2_i].trg_n & 0x2)
                            push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[f2_i].pre_trg);
                    } else if (f_msg->fa_msg[f1_i].trg_n & 0x2)
                        push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[f1_i].pre_trg);
                    // merged into "res_split"
                    //cut_cigar(split_cigar, split_len, la, APP->read_len);
                    la->split_flag = 1;
                } else if (res & 2) la->split_flag = 1;
                //else {printf("del-score: 0\t"); printcigar(stdout, *split_cigar, *split_len); printf("\n");}
            }
            //else	//dis > split_read_len
        } else if (dis < -APP->match_dis) {
#ifdef __DEBUG__
            fprintf(stdout, "%d\t%lld\t%d\tINS\n", a_msg[s1_i].at[s1_aln_i].chr, (long long)pos, (0-dis));
#endif
			split_ref_len = split_read_len + dis;
            if (split_ref_len < 0) // overlap ins
            {
                fprintf(stderr, "\n[split_map] error breakpoint.\n");
                exit(0);
            }
            else 
            {
                split_ref_seq = (uint8_t*)malloc(split_ref_len * sizeof(uint8_t));
                ref_offset = a_msg[s1_i].at[s1_aln_i].offset + APP->seed_len + a_msg[s1_i].at[s1_aln_i].len_dif;
                pac2fa_core(bns, pac, a_msg[s1_i].at[s1_aln_i].chr, ref_offset-1, &split_ref_len, 1, &N_flag, &N_len, split_ref_seq);
                if (split_ref_len < hash_len) 
                {
                    cigar32_t *k_cigar=0;
                    int k_clen;
                    int score = ksw_global(split_read_len, split_read_seq, split_ref_len, split_ref_seq, 5, bwasw_sc_mat, 5, 2, abs(split_ref_len-split_read_len)+3, &k_clen, &k_cigar);
                    //printf("score3: %d\n", score);
                    *split_len = 0;
                    _push_cigar(split_cigar, split_len, split_m, k_cigar, k_clen);
                    free(k_cigar);
                }
                else
                {
                    res = split_indel_map(split_cigar, split_len, split_m, split_read_seq, split_read_len, split_ref_seq, split_ref_len, ref_offset, hash_len, hash_step, AP->split_len, hash_num, hash_node, key_len, hash_size);
                    if (res & 1)
                    {
                        if (nsrand == 1) {
                            if (f_msg->fa_msg[f2_i].trg_n & 0x2)
                                push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[f2_i].pre_trg);
                        } else if (f_msg->fa_msg[f1_i].trg_n & 0x2)
                                push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[f1_i].pre_trg);
                           
                        //cut_cigar(split_cigar, split_len, la, APP->read_len);
                        la->split_flag = 1;
                    } else if (res & 2) la->split_flag = 1;
                }
            }
		}
		//for MIS-MATCH
		else {
			//XXX	MIS-MATCH
			split_ref_len = split_read_len + dis;
			split_ref_seq = (uint8_t*)malloc(split_ref_len * sizeof(uint8_t));
			ref_offset = a_msg[s1_i].at[s1_aln_i].offset + APP->seed_len + a_msg[s1_i].at[s1_aln_i].len_dif;
			pac2fa_core(bns, pac, a_msg[s1_i].at[s1_aln_i].chr, ref_offset-1, &split_ref_len, 1, &N_flag, &N_len, split_ref_seq);
			int min_b = abs(dis)+3;
			cigar32_t *cigar=0; int c_len;
			int score = ksw_global(split_read_len, split_read_seq, split_ref_len, split_ref_seq, 5, bwasw_sc_mat, 5, 2, min_b, &c_len, &cigar);

			//XXX 'spurious MATCH'
			if (score < 0) 
			{
                res = hash_both_bound_map(split_cigar, split_len, split_m, split_ref_seq, split_ref_len, split_read_seq, split_read_len, 
                                        hash_num, hash_node, AP->hash_len, AP->hash_key_len, AP->hash_step, AP->split_len);
                if (res & 1)//XXX pull trigger, multi triggers
				{
					if (nsrand == 1) {
                        if (f_msg->fa_msg[f2_i].trg_n & 0x2)
                            push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[f2_i].pre_trg);
                    } else if (f_msg->fa_msg[f1_i].trg_n & 0x2)
                        push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[f1_i].pre_trg);
                    //cut aln-cigar into split-cigar on 'nSmH' 
                    //cut_cigar(split_cigar, split_len, la, APP->read_len);
                    la->split_flag = 1;
                } else if (res & 2) la->split_flag = 1;
            }
			else
			{
                *split_len = 0;
                _push_cigar(split_cigar, split_len, split_m, cigar, c_len);
			}
			free(cigar);
		}
	}
	else	//dif srand
	{
        int pos1, pos2;
        if (a_msg[s1_i].at[s1_aln_i].nsrand == 1)
            pos1 = a_msg[s1_i].at[s1_aln_i].offset+APP->seed_len-1;
        else pos1 = a_msg[s1_i].at[s1_aln_i].offset;
        if (a_msg[s2_i].at[s2_aln_i].nsrand == 1)
            pos2 = a_msg[s2_i].at[s2_aln_i].offset;
        else pos2 = a_msg[s2_i].at[s2_aln_i].offset+APP->seed_len-1;
		fprintf(stderr, "\n%d\t%d\tcomplex SV\n", a_msg[s1_i].at[s1_aln_i].chr, pos1);
		fprintf(stderr, "%d\t%d\tcomplex SV\n", a_msg[s2_i].at[s2_aln_i].chr, pos2);
	}
    _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), *split_cigar, *split_len);
	free(split_read_seq); free(split_ref_seq);
}

void check_cigar(cigar32_t *cigar, int cigar_len, char *read_name, int read_len)
{
	int i;
	int mid[5] = {0, 0, 0, 0, 0};
	for (i = 0; i < cigar_len; ++i)
		mid[(int)(cigar[i]&0xf)] += (cigar[i]>>4);
	fprintf(stdout, "\t%d M, %d I, %d D, %d S\t", mid[0], mid[1], mid[2], mid[4]);
    if (mid[0] + mid[1] + mid[4] != read_len) { fprintf(stderr, "\n[cigar len error]: %s cigar error[%d]\t", read_name, read_len); printcigar(stderr, cigar, cigar_len); }
}

//do hash map	XXX
int frag_head_bound_fix(frag_msg *f_msg, aln_msg *a_msg, 
						cigar32_t **cigar, int *cigar_len, int *cigar_m,
						uint64_t *offset, bntseq_t *bns, uint8_t *pac, 
						char *read_seq, uint8_t *seq1, uint8_t *seq2, 
                        lsat_aln_per_para *APP, lsat_aln_para *AP,
						uint32_t **hash_num, uint64_t ***hash_node, 
						line_aln_res *la)//aln_res *a_res)
{
	int left_bound = f_msg->frag_left_bound;
	int frag_i, seed_x, seed_i, aln_i, i;
	int N_flag, N_len;
	int read_len, ref_len, read_start/*0-base*/;
	uint64_t ref_start;
	if (f_msg->fa_msg[0].srand == 1) {	//'+' srand
		frag_i = f_msg->frag_num-1;
		seed_x = f_msg->fa_msg[frag_i].seed_num - 1;
		seed_i = f_msg->fa_msg[frag_i].seed_i[seed_x];
		aln_i = f_msg->fa_msg[frag_i].seed_aln_i[seed_x];
		if (a_msg[seed_i].read_id != 1) {
		//if ((a_msg[seed_i].read_id - left_bound) > 1)
			read_len = ((left_bound==0)?0:APP->seed_inv) + (a_msg[seed_i].read_id - left_bound - 1) * APP->seed_step;
            if (read_len < 0) { 
                fprintf(stderr, "[frag_head_bound] Error: %s left bound error.\n", APP->read_name);     
                exit(0); 
            }
			read_start = (left_bound==0)?0:(left_bound * APP->seed_step - APP->seed_inv);
		} else {
			(*cigar_len) = 0;
			*offset = a_msg[seed_i].at[aln_i].offset;
            la->res[la->cur_res_n].offset = a_msg[seed_i].at[aln_i].offset;
            la->res[la->cur_res_n].cigar_len = 0;
			return 0;
		}
	}
	else {	//'-' srand
		frag_i = 0; seed_x = 0;
		seed_i = f_msg->fa_msg[0].seed_i[0];
		aln_i = f_msg->fa_msg[0].seed_aln_i[0];
		read_len = ((left_bound==0)?APP->last_len:APP->seed_inv) + (a_msg[seed_i].read_id - 1 - left_bound) * APP->seed_step;
        if (read_len < 0) { 
            fprintf(stderr, "[frag_head_bound] Error: %s left bound error.\n", APP->read_name); 
            exit(0); 
        }
		read_start= (left_bound==0)?0:(APP->last_len+left_bound*APP->seed_step-APP->seed_inv);
	}
	//hash map
	{	
		//XXX un-overlap
		int hash_len = AP->hash_len, hash_step = AP->hash_step, hash_key = AP->hash_key_len;
		(*offset) = a_msg[seed_i].at[aln_i].offset;
        la->res[la->cur_res_n].offset = a_msg[seed_i].at[aln_i].offset;
		ref_len = read_len + hash_step * 2;	//XXX
		for (i = 0; i < read_len; ++i)
			seq1[i] = nst_nt4_table[(int)read_seq[read_start+i]];
		//ref
		ref_start = a_msg[seed_i].at[aln_i].offset - ref_len;	//1-base
		//ref XXX check for 'N'
		pac2fa_core(bns, pac, a_msg[seed_i].at[aln_i].chr, ref_start-1/*0-base*/, &ref_len, 1, &N_flag, &N_len, seq2);

        line_node trg;
        int res;
        res = hash_left_bound_map(cigar, cigar_len, cigar_m, seq2, ref_len, seq1, read_len, hash_num, hash_node, hash_len, hash_key, hash_step, AP->split_len);
		if (res & 1) {
			if (f_msg->fa_msg[0].srand == 1) { 	//'+' srand
                if (f_msg->fa_msg[frag_i].trg_n & 0x2) 
                    push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[frag_i].pre_trg);
            } else if (f_msg->fa_msg[0].trg_n & 0x1)
                push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[0].next_trg);
            //XXX
            la->split_flag = 1;
		} else if (res & 2) la->split_flag = 1;
		for (i = 0; i < *cigar_len; ++i) {
			if (((*cigar)[i] & 0xf) == CMATCH || ((*cigar)[i] & 0xf) == CDEL) {
				(*offset) -= ((*cigar)[i]>>4); 
                (la->res[la->cur_res_n].offset) -= ((*cigar)[i]>>4);
            }
        }
        _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), *cigar, *cigar_len);
	}
	return 0;
}

int frag_tail_bound_fix(frag_msg *f_msg, aln_msg *a_msg, 
						cigar32_t **cigar, int *cigar_len, int *cigar_m,
						bntseq_t *bns, uint8_t *pac, 
						char *read_seq, uint8_t *seq1, uint8_t *seq2, 
                        lsat_aln_per_para *APP, lsat_aln_para *AP,
						uint32_t **hash_num, uint64_t ***hash_node, 
						line_aln_res *la)
{
	int right_bound = f_msg->frag_right_bound;
	int frag_i, seed_x, seed_i, aln_i, i;
	int N_flag, N_len;
	int read_len, read_start, ref_len;
	uint64_t ref_start;

	if (f_msg->fa_msg[0].srand == 1) {	//'+' srand
		seed_i = f_msg->fa_msg[0].seed_i[0];
		aln_i = f_msg->fa_msg[0].seed_aln_i[0];
		read_start = a_msg[seed_i].read_id * APP->seed_step - APP->seed_inv;
		read_len = ((right_bound==APP->seed_all+1)?APP->last_len:APP->seed_inv) + (right_bound - 1 - a_msg[seed_i].read_id) * APP->seed_step;
        if (read_len < 0) { fprintf(stderr, "[frag_tail_bound] Error: %s right bound error.\n", APP->read_name); exit(0); }
	} else {
		frag_i = f_msg->frag_num - 1;
		seed_x = f_msg->fa_msg[frag_i].seed_num - 1;
		seed_i = f_msg->fa_msg[frag_i].seed_i[seed_x];
		aln_i = f_msg->fa_msg[frag_i].seed_aln_i[seed_x];
		if (a_msg[seed_i].read_id != APP->seed_all)
		{
			read_start = a_msg[seed_i].read_id * APP->seed_step - APP->seed_inv + APP->last_len;
			read_len = ((right_bound== APP->seed_all+1)?0:APP->seed_inv) + (right_bound - 1 - a_msg[seed_i].read_id) * APP->seed_step;
            if (read_len < 0) { fprintf(stderr, "[frag_tail_bound] Error: %s right bound error.\n", APP->read_name); exit(0); }
		}
		else
		{
			(*cigar_len) = 0;
			return 0;
		}
	}
	//hash map
	{
		//XXX un-overlap
		int hash_len = AP->hash_len, hash_step = AP->hash_step, hash_key = AP->hash_key_len;
		ref_len = read_len + hash_step * 2;
		//read
		for (i = 0; i < read_len; ++i)
			seq1[i] = nst_nt4_table[(int)read_seq[read_start+i]];
		//ref
		ref_start = a_msg[seed_i].at[aln_i].offset + APP->seed_len+ a_msg[seed_i].at[aln_i].len_dif;	//1-base
		//printf("tail, ref-start: %lld, read-start: %d\n", ref_start, read_start);
		pac2fa_core(bns, pac, a_msg[seed_i].at[aln_i].chr, ref_start-1/*0-base*/, &ref_len, 1, &N_flag, &N_len, seq2);

		//fprintf(stderr, "read: %s\nref: %s\n", read, ref);
        line_node trg;
        int res;
        res = hash_right_bound_map(cigar, cigar_len, cigar_m, seq2, ref_len, seq1, read_len, hash_num, hash_node, hash_len, hash_key, hash_step, AP->split_len);
        if (res & 1)	//XXX pull trigger
		{
            // pull trigger
			if (f_msg->fa_msg[0].srand == 1) { //'+' srand
                if (f_msg->fa_msg[0].trg_n & 0x1)
                    push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[0].next_trg);
            } else if (f_msg->fa_msg[frag_i].trg_n & 0x2)
                push_trg(&(la->trg), &(la->trg_n), &(la->trg_m), f_msg->fa_msg[frag_i].pre_trg);
            //cut aln-cigar into split-cigar on 'nSmH' 
            //cut_cigar(cigar, cigar_len, la, APP->read_len); 
            la->split_flag = 1;
		} else if (res & 2) la->split_flag = 1;
        _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), *cigar, *cigar_len);
	}
	return 0;
}

// XXX need a split-flag, generated by upstream program?
// generate split alignment, if necessary
    // 1. fix split_i : tail_S
    // 2. create split_i+1: offset, cigar, head_S, push_cigar
void lsat_res_split(line_aln_res *la, int read_len, lsat_aln_para *AP)
{
    int i, j, res_n, s_i;
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
            if ((len = cigar[j] >> 4) >= AP->split_len) { // ins split
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
            if ((len = cigar[j] >> 4) >= AP->split_len) { // del split
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
                // tail S
                len1 = readInCigar(la->res[res_n].cigar, la->res[res_n].cigar_len);
                tail_s = read_len - len1;
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((tail_s << 4) | CSOFT_CLIP));
                // ++res_n
                ++res_n; push_res(la);
                // offset
                la->res[res_n].offset = la->res[res_n-1].offset + refInCigar(la->res[res_n-1].cigar, la->res[res_n-1].cigar_len) + (int)(cigar[j+1] >> 4);
                // head S
                len = cigar[j] >> 4;
                head_s = len1 + len;
                _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), (cigar32_t)((head_s << 4) | CSOFT_CLIP));
            } else _push_cigar1(&(la->res[res_n].cigar), &(la->res[res_n].cigar_len), &(la->res[res_n].c_m), cigar[j]);
        } else if (op != CHARD_CLIP) {
            fprintf(stderr, "[lsat_res_split] Error: Unexpected cigar operation: %d.\n", op);
            exit(0);
        }
    }
    free(cigar);
}
// generate NM/AS tags
void lsat_res_aux(line_aln_res *la, bntseq_t *bns, uint8_t *pac, char *read_seq, int read_len, lsat_aln_para *AP, lsat_aln_per_para *APP)
{
    //int read_len = a_res->read_len;
    char *ref_seq=0; int ref_len, N_flag, N_len;
    int m, i, j, k, op, len, ref_i, read_i;
    int n_mm, n_m, n_o, n_e; // mis-match/match/gap-open/gap-extension
    int mm_tmp;
    res_t *r; cigar32_t *cigar;

    //for (l = 0; l < a_res->l_m; ++l) {
    for (m = 0; m <= la->cur_res_n; ++m) {
        r = la->res+m; 
        cigar = r->cigar;
        ref_len = refInCigar(cigar, r->cigar_len);
        ref_seq = (char*)realloc(ref_seq, ref_len * sizeof(char));
        pac2fa(bns, pac, r->chr, r->offset-1/*0-base*/, &ref_len, 1, &N_flag, &N_len, ref_seq);
        ref_i = read_i = n_mm = n_m = n_o = n_e = 0;
        for (i = 0; i < r->cigar_len; ++i) {
            op = cigar[i] & 0xf, len = cigar[i]>>4;
            if (op == CMATCH) {
                mm_tmp = 0;
                for (j = 0; j < len; ++j) {
                    if (read_seq[read_i++] != ref_seq[ref_i++])
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
            } else { fprintf(stderr, "[lsat_gen_aux] Error: Unexpected cigar operation: %d.\n", op); 
                printcigar(stderr, cigar, r->cigar_len); 
                exit(0); }
        }
        if (read_i != read_len || ref_i != ref_len) { fprintf(stderr, "[lsat_gen_aux] Error: %s Unmatched length: read: %d %d.\tref: %d %d\n", APP->read_name, read_i, read_len, ref_i, ref_len);
            printcigar(stderr, cigar, r->cigar_len); fprintf(stderr, "\n");
            exit(0); }
        // calculate NM and AS
        r->NM = n_mm + n_e;
        r->score = n_m * AP->match - n_mm * AP->mis - n_o * AP->gapo - n_e * AP->gape;
        la->tol_score += r->score;
        la->tol_NM += r->NM;
    }
    la->tol_score -= (la->cur_res_n * AP->split_pen);
    //}
    free(ref_seq);
}
// filter OPT-res or multi-OPT, based on merged-message
    // before filtering:
    // (L_NMERG, x) // not merged
    // (L_MERGH, x) // Merged, head
    // (L_MERGB, i) // Merged, body, head=i
    // (L_INTER, x) // Inter line
        //  0,x -> Merged and Head -> b,s
        // -1,x -> Not
        // -2,i -> Merged, Head is i
    // after filtering:
    //  0/-1/-2,x: dump OR is secondary
    //  1,s: keep and is best, secondary is [s], -1 for NULL
void res_filter(aln_res *a_res)
{
    int i, l, m_i, m_ii;  // number of merged block 
    line_node *m_b; // merged line #, {x,y}, x:best, y:secondary, (-1 for NULL)
    int *head;

    l = a_res->l_m; m_i = -1;
    m_b = (line_node*)malloc(l * sizeof(line_node));
    head = (int*)malloc(l * sizeof(line_node));
    for (i = 0; i < l; ++i) {
        if (a_res->la[i].merg_msg.x & L_NMERG)
            a_res->la[i].merg_msg = (line_node){1,-1};
        else if (a_res->la[i].merg_msg.x & L_MERGH) {
            ++m_i;
            m_b[m_i].x = i, m_b[m_i].y = -1; // init of m_bh
            head[m_i] = i;
        } else { // Merged, check for Best and Secondary in m_b[m_i];
            // get head: m_ii
            for (m_ii = 0; m_ii <= m_i; ++m_ii) {
                if (head[m_ii] == a_res->la[i].merg_msg.y) break;
            }
            if (m_ii > m_i) { 
                fprintf(stderr, "[resust_filter] %s BUG.\n", READ_NAME); exit(-1); 
            }
            if (a_res->la[i].tol_score > a_res->la[m_b[m_ii].x].tol_score
            || (a_res->la[i].tol_score == a_res->la[m_b[m_ii].x].tol_score 
             && a_res->la[i].tol_NM < a_res->la[m_b[m_ii].x].tol_NM)) { // new best
                m_b[m_ii].y = m_b[m_ii].x; m_b[m_ii].x = i;
            } 
            else if (m_b[m_ii].y == -1) m_b[m_ii].y = i;
            else if (a_res->la[i].tol_score > a_res->la[m_b[m_ii].y].tol_score
              || (a_res->la[i].tol_score == a_res->la[m_b[m_ii].y].tol_score
               && a_res->la[i].tol_NM < a_res->la[m_b[m_ii].y].tol_NM)) {
              m_b[m_ii].y = i;
            } 
        }
    }
    // calculate mapping quality by best and secondary aln XXX
    for (i = 0; i <= m_i; ++i) {
        if (m_b[i].y >= 0 && a_res->la[m_b[i].y].tol_score > (a_res->la[m_b[i].x].tol_score) / 2) {//XXX
            a_res->la[m_b[i].x].merg_msg = (line_node){1,m_b[i].y};
        } else a_res->la[m_b[i].x].merg_msg = (line_node){1,-1};
    }
    free(m_b); free(head);
}

//read_seq: char or uint8_t?
void frag_check(aln_msg *a_msg, frag_msg **f_msg, aln_res *a_res,
                    bntseq_t *bns, uint8_t *pac, 
                    const char *read_prefix, char *read_seq,
                    lsat_aln_per_para *APP, lsat_aln_para *AP,
                    int line_n, uint32_t **hash_num, uint64_t ***hash_node)  //XXX seed_len -> seed_len+inv :OK
{
    // for debug
    strcpy(READ_NAME, APP->read_name);
	int i, j;
	int max_len = APP->read_len;
	uint8_t *seq1 = (uint8_t*)malloc((max_len+1)*sizeof(uint8_t));
	uint8_t *seq2 = (uint8_t*)malloc((max_len+1)*sizeof(uint8_t));
	
    // realloc mem, if necessary
    if (line_n > a_res->l_m) {
        if ((a_res->la = (line_aln_res*)realloc(a_res->la, line_n * sizeof(line_aln_res))) == NULL) {
            fprintf(stderr, "[frag_check] Not enough memory.\n"); exit(0);
        }
        for (i = a_res->l_m; i < line_n; ++i) {
            a_res->la[i].res_m = 10, a_res->la[i].cur_res_n = 0, a_res->la[i].split_flag = 0;
            a_res->la[i].res = (res_t*)malloc(10 * sizeof(res_t));
            for (j = 0; j < 10; ++j) {
                a_res->la[i].res[j].c_m = 100;
                a_res->la[i].res[j].cigar = (cigar32_t*)malloc(100 * sizeof(cigar32_t));
                a_res->la[i].res[j].cigar_len = 0;
            }
            a_res->la[i].tol_score = a_res->la[i].tol_NM = 0;
            a_res->la[i].trg_m = 10, a_res->la[i].trg_n = 0;
            a_res->la[i].trg = (line_node*)malloc(10 * sizeof(line_node));
        }
        a_res->l_m = line_n;
    }
    // initialization
    a_res->l_n = line_n;
    for (i = 0; i < line_n; ++i) {
        a_res->la[i].merg_msg = ((*f_msg)+i)->merg_msg;
        a_res->la[i].cur_res_n = 0;
        a_res->la[i].split_flag = 0;
        a_res->la[i].tol_score = 0; a_res->la[i].tol_NM = 0;
        a_res->la[i].trg_n = 0;
        for (j = 0; j < a_res->la[i].res_m; ++j) a_res->la[i].res[j].cigar_len = 0;
    }

    cigar32_t *cigar;
    int cigar_len, cigar_m, res_n;
    uint64_t offset;
    cigar_m = 100;
    cigar_len = 0;
    cigar = (cigar32_t*)malloc(cigar_m * sizeof(cigar32_t));
    // check for every line
    for (j = 0; j < line_n; ++j) {
        a_res->la[j].res[a_res->la[j].cur_res_n].cigar_len = 0;
        a_res->la[j].res[a_res->la[j].cur_res_n].nsrand = (((*f_msg)+j)->fa_msg[0].srand == 1)?1:0;
        a_res->la[j].res[a_res->la[j].cur_res_n].chr = ((*f_msg)+j)->fa_msg[0].chr;
		/* extend once
		   for every frag : take it as a RIGHT frag
		   for every interval : SV breakpoint(s).*/
		if (((*f_msg)+j)->fa_msg[0].srand == 1) {
			//head 'S'
				if (((*f_msg)+j)->frag_left_bound > 0) {
					cigar_len = 1;
					cigar[0] = ((((*f_msg)+j)->frag_left_bound * APP->seed_step - APP->seed_inv) << 4) | CSOFT_CLIP;
                    res_n = a_res->la[j].cur_res_n;
                    _push_cigar(&(a_res->la[j].res[res_n].cigar), &(a_res->la[j].res[res_n].cigar_len), &(a_res->la[j].res[res_n].c_m), cigar, cigar_len);
				}
			//fix the boundary blank before first frag, if it exists
                frag_head_bound_fix((*f_msg)+j, a_msg, &cigar, &cigar_len, &cigar_m, &offset, bns, pac, read_seq, seq1, seq2, APP, AP, hash_num, hash_node, &(a_res->la[j]));
                for (i = ((*f_msg)+j)->frag_num-1; i > 0; --i) {
                    frag_extend((*f_msg)+j, a_msg, i, bns, pac, read_seq, APP->read_len, seq1, seq2, APP, AP, &(a_res->la[j]));
                    split_mapping(&cigar, &cigar_len, &cigar_m, bns, pac, read_seq, ((*f_msg)+j), a_msg, hash_num, hash_node, i, i-1, APP, AP, &(a_res->la[j]));
                }
                frag_extend((*f_msg)+j, a_msg, i, bns, pac, read_seq, APP->read_len, seq1, seq2, APP, AP, &(a_res->la[j]));
			//fix the boundary blank after the last frag
                frag_tail_bound_fix((*f_msg)+j, a_msg, &cigar, &cigar_len, &cigar_m, bns, pac, read_seq, seq1, seq2, APP, AP, hash_num, hash_node, &(a_res->la[j]));
			//tail 'S'
				if (((*f_msg)+j)->frag_right_bound <= APP->seed_all) {
					cigar_len = 1;
                    cigar[0] = ((APP->read_len - (((*f_msg)+j)->frag_right_bound-1) * APP->seed_step) << 4) | CSOFT_CLIP; 
                    res_n = a_res->la[j].cur_res_n;
                    _push_cigar(&(a_res->la[j].res[res_n].cigar), &(a_res->la[j].res[res_n].cigar_len), &(a_res->la[j].res[res_n].c_m), cigar, cigar_len);
				}
            //split and aux
            lsat_res_split(a_res->la+j, APP->read_len, AP);
            lsat_res_aux(a_res->la+j, bns, pac, read_seq, APP->read_len, AP, APP);
		}
		else {	//'-' srand
			//convert into rev-com
                // XXX convert ONLY once during all processes?
                char *reco_read_seq = (char*)malloc((APP->read_len+1) * sizeof(char));
                for (i = 0; i < APP->read_len; ++i) 
                    reco_read_seq[i] = (read_seq[APP->read_len-1-i]=='A')?'T':((read_seq[APP->read_len-1-i]=='C')?'G':
                            ((read_seq[APP->read_len-1-i]=='G')?'C':((read_seq[APP->read_len-1-i]=='T')?'A':'N')));
                reco_read_seq[i] = 0;
                for (i = 0; i < APP->seed_all; ++i) a_msg[i].read_id = (APP->seed_all + 1 - a_msg[i].read_id);
                int tmp = ((*f_msg)+j)->frag_left_bound;
                ((*f_msg)+j)->frag_left_bound = APP->seed_all + 1 - ((*f_msg)+j)->frag_right_bound;
                ((*f_msg)+j)->frag_right_bound = APP->seed_all + 1 - tmp;
			//head 'S'
				if (((*f_msg)+j)->frag_left_bound > 0) {
					cigar_len = 1;
					cigar[0] = ((((*f_msg)+j)->frag_left_bound * APP->seed_step - APP->seed_inv + APP->last_len) << 4) | CSOFT_CLIP;
                    res_n = a_res->la[j].cur_res_n;
                    _push_cigar(&(a_res->la[j].res[res_n].cigar), &(a_res->la[j].res[res_n].cigar_len), &(a_res->la[j].res[res_n].c_m), cigar, cigar_len);
				}
			//fix the boundary blank before first frag
                frag_head_bound_fix((*f_msg)+j, a_msg, &cigar, &cigar_len, &cigar_m, &offset, bns, pac, reco_read_seq, seq1, seq2, APP, AP, hash_num, hash_node, &(a_res->la[j]));
                for (i = 0; i < ((*f_msg)+j)->frag_num-1; ++i) {
                    frag_extend((*f_msg)+j, a_msg, i, bns, pac, reco_read_seq, APP->read_len, seq1, seq2, APP, AP, &(a_res->la[j]));
                    split_mapping(&cigar, &cigar_len, &cigar_m, bns, pac, reco_read_seq, (*f_msg)+j, a_msg, hash_num, hash_node, i, i+1, APP, AP, &(a_res->la[j]));
                }
                frag_extend(((*f_msg)+j), a_msg, i, bns, pac, reco_read_seq, APP->read_len, seq1, seq2, APP, AP, &(a_res->la[j]));
			//fix the boundary blank after the last frag
                frag_tail_bound_fix(((*f_msg)+j), a_msg, &cigar, &cigar_len, &cigar_m, bns, pac, reco_read_seq, seq1, seq2, APP, AP, hash_num, hash_node, &(a_res->la[j]));
			//tail 'S'
				if (((*f_msg)+j)->frag_right_bound <= APP->seed_all) {
					cigar_len = 1;
                    cigar[0] = (((APP->seed_all - ((*f_msg)+j)->frag_right_bound+1) * APP->seed_step - APP->seed_inv) << 4) | CSOFT_CLIP;
                    res_n = a_res->la[j].cur_res_n;
                    _push_cigar(&(a_res->la[j].res[res_n].cigar), &(a_res->la[j].res[res_n].cigar_len), &(a_res->la[j].res[res_n].c_m), cigar, cigar_len);
				}
            lsat_res_split(a_res->la+j, APP->read_len, AP);
            lsat_res_aux(a_res->la+j, bns, pac, reco_read_seq, APP->read_len, AP, APP);
			free(reco_read_seq);
			for (i = 0; i < APP->seed_all; ++i) a_msg[i].read_id = (APP->seed_all + 1 - a_msg[i].read_id);
		}
	}
    // filter line-aln-res in merged-lines, get opt-res OR multi-sub-opt-res
    res_filter(a_res);

	free(seq1); free(seq2); free(cigar);
}
