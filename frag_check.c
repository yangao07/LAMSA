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

void frag_init_msg(frag_msg *f_msg, int frag_max)
{
	f_msg->frag_max = frag_max; 
	f_msg->frag_num = 0;
	f_msg->per_seed_max = frag_max;
	f_msg->seed_num = 0;

	f_msg->fa_msg = (frag_aln_msg*)malloc(frag_max*sizeof(frag_aln_msg));
	int i;
	for (i = 0; i < frag_max; ++i)
	{
		f_msg->fa_msg[i].cigar_max = 100;
		f_msg->fa_msg[i].cigar = (uint32_t*)malloc(100 * sizeof(uint32_t));
		f_msg->fa_msg[i].cigar_len = 0;

		f_msg->fa_msg[i].seed_i = (int*)malloc(f_msg->per_seed_max*sizeof(int));
		f_msg->fa_msg[i].seed_aln_i = (int*)malloc(f_msg->per_seed_max*sizeof(int));
		f_msg->fa_msg[i].seed_num = 0;
        f_msg->fa_msg[i].next_trigger_n = 0;
        f_msg->fa_msg[i].pre_trigger_n = 0;
        f_msg->fa_msg[i].trigger_m = 100;
        f_msg->fa_msg[i].trigger = (int*)malloc(100 * sizeof(int));
	}
}

int frag_copy_msg(frag_msg *ff_msg, frag_msg *tf_msg)
{
	tf_msg->frag_max = ff_msg->frag_max;
	tf_msg->frag_num = 0;
	tf_msg->last_len = ff_msg->last_len;
	tf_msg->seed_num = 0;
	tf_msg->seed_all = ff_msg->seed_all;
	tf_msg->per_seed_max = ff_msg->per_seed_max;
	//tf_msg->fa_msg = (frag_aln_msg*)malloc(tf_msg->frag_max * sizeof(frag_aln_msg));
	int i;
	for (i = 0; i < tf_msg->frag_max; ++i)
	{
		tf_msg->fa_msg[i].cigar_max = 100;
		//tf_msg->fa_msg[i].cigar = (uint32_t*)malloc(100 * sizeof(uint32_t));
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
			free(f_msg[i].fa_msg[j].trigger);
		}
		free(f_msg[i].fa_msg);
	}
	free(f_msg);
}

int frag_set_msg(aln_msg *a_msg, int seed_i, int aln_i, 
				 int FLAG, 
				 frag_msg *f_msg, int frag_i, 
				 int seed_len)//FLAG 0:start / 1:end / 2:seed
{
	if (frag_i == 0)
		f_msg->seed_num=0;
	if (FLAG == FRAG_END)	//end
	{
		f_msg->fa_msg[frag_i].chr = a_msg[seed_i].at[aln_i].chr;
		f_msg->fa_msg[frag_i].srand = a_msg[seed_i].at[aln_i].nsrand;	//+:1/-:-1
		f_msg->fa_msg[frag_i].seed_i[0] = seed_i;
		f_msg->fa_msg[frag_i].seed_aln_i[0] = aln_i;
		f_msg->fa_msg[frag_i].seed_num = 1;
		f_msg->seed_num++;
		f_msg->fa_msg[frag_i].flag = UNCOVERED;
		f_msg->fa_msg[frag_i].cigar_len = 0;
	}
	else if(FLAG==FRAG_START)	//start
	{
		f_msg->frag_num = frag_i + 1;
	}
	else			//seed
	{
		f_msg->fa_msg[frag_i].seed_i[f_msg->fa_msg[frag_i].seed_num] = seed_i;
		f_msg->fa_msg[frag_i].seed_aln_i[f_msg->fa_msg[frag_i].seed_num] = aln_i;
		f_msg->seed_num++;
		f_msg->fa_msg[frag_i].seed_num++;
	}
	return 0;
}

int frag_trigger_set(frag_dp_node f_node, frag_msg *f_msg, int frag_i)
{
    f_msg->fa_msg[frag_i].next_trigger_n = f_node.next_trigger_n;
    f_msg->fa_msg[frag_i].pre_trigger_n = f_node.pre_trigger_n;
    int i;
    if (f_node.next_trigger_n + f_node.pre_trigger_n > f_msg->fa_msg[frag_i].trigger_m)
    {
        f_msg->fa_msg[frag_i].trigger = (int*)realloc(f_msg->fa_msg[frag_i].trigger, (f_node.next_trigger_n+f_node.pre_trigger_n) * sizeof(int));
        f_msg->fa_msg[frag_i].trigger_m = f_node.next_trigger_n+f_node.pre_trigger_n;
    }
    for (i = 0; i < f_node.next_trigger_n+f_node.pre_trigger_n; ++i)
        f_msg->fa_msg[frag_i].trigger[i] = f_node.trigger[i];
    return 0;
}

/*int frag_refresh(frag_msg* f_msg, int f_i, int ref_offset, int read_offset, int FLAG)
{
	if (FLAG == RIGHT)	//right offset
	{
		f_msg->fa_msg[f_i].ex_read_end += read_offset;
		if (f_msg->fa_msg[f_i].srand == 1)	//+ srand
			f_msg->fa_msg[f_i].ex_ref_end += ref_offset;
		else						//- srand
			f_msg->fa_msg[f_i].ex_ref_end -= ref_offset;
		//f_msg->fa_msg[f_i].per_n += read_offset / PER_LEN;
	}
	else				//left offset
	{
		f_msg->fa_msg[f_i].ex_read_begin -= read_offset;
		if (f_msg->fa_msg[f_i].srand == 1)	//+ srand
			f_msg->fa_msg[f_i].ex_ref_begin -= ref_offset;
		else
			f_msg->fa_msg[f_i].ex_ref_begin += ref_offset;
		//f_msg->fa_msg[f_i].per_n += read_offset / PER_LEN;
	}
	//printf("Refresh: %d %d %d %d\n", f_i, ref_offset, read_offset, FLAG);
	return 0;
}*/

/*int frag_covered(frag_msg *f_msg, int f_i, int b_f)
{
	f_msg->fa_msg[b_f].flag = COVERED;
	f_msg->fa_msg[f_i].b_f = b_f;
	//	printf("Covered: %d C %d\n", f_i, b_f);
	return 0;
}*/

/*int unsoap2dp_extend(frag_msg *f_msg, int f_i, int FLAG, bntseq_t *bns, uint8_t *pac, char *read_seq, uint8_t *seq1, uint8_t *seq2)
{
	int ret;
	int ref_len, read_len;
	int ref_l_os, read_l_os, ref_r_os, read_r_os;
	int b_f = f_msg->fa_msg[f_i].b_f;
	int ff, m;
	char *seq_s2 = malloc(16385*sizeof(char));

	// XXX Is it necessary to take a 'COVERED frag' as an individaul frag in the next reverse extend?
	//	   Here is 'IS'.
	
	if (FLAG == RIGHT)
	{
		if (b_f == 0)
		{
			read_len = f_msg->last_len; //XXX
			ref_len = read_len;
		}
		else
		{
			read_len = f_msg->fa_msg[b_f-1].read_begin - f_msg->fa_msg[b_f].read_end - 1;
			ref_len = read_len;
		}
		//pac2fa: for '-' srand, return the re-co seq(same orientation as read)
		if (f_msg->fa_msg[b_f].srand == -1)	//- srand
			pac2fa_core(bns, pac, f_msg->fa_msg[b_f].chr, f_msg->fa_msg[b_f].ref_end-1-ref_len, &ref_len, f_msg->fa_msg[b_f].srand, &ff, seq1);
		else
			pac2fa_core(bns, pac, f_msg->fa_msg[b_f].chr, f_msg->fa_msg[b_f].ref_end, &ref_len, f_msg->fa_msg[b_f].srand, &ff, seq1);
		strncpy(seq_s2, read_seq+f_msg->fa_msg[b_f].read_end, read_len);
		for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
		ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
		frag_refresh(f_msg, f_i, ref_r_os, read_r_os, RIGHT);
	}
	else	//LEFT
	{
		if (b_f == f_msg->frag_num - 1)
		{
			if (f_msg->fa_msg[b_f].read_begin != 1)	//Un-soap2dp HEAD
			{
				read_len = f_msg->fa_msg[b_f].read_begin - 1;
				ref_len =read_len;
			}
			else
				return BAD_ALIGN;
		}
		else
		{
			read_len = f_msg->fa_msg[b_f].read_begin - f_msg->fa_msg[b_f+1].read_end - 1;
			ref_len = read_len;
		}
		if (f_msg->fa_msg[b_f].srand == -1)
			pac2fa_core(bns, pac, f_msg->fa_msg[b_f].chr, f_msg->fa_msg[b_f].ref_begin, &ref_len, f_msg->fa_msg[b_f].srand, &ff, seq1);
		else
			pac2fa_core(bns, pac, f_msg->fa_msg[b_f].chr, f_msg->fa_msg[b_f].ref_begin-ref_len-1, &ref_len, f_msg->fa_msg[b_f].srand, &ff, seq1);
		strncpy(seq_s2, read_seq+(f_msg->fa_msg[b_f].read_begin - read_len -1), read_len);
		for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
		ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
		frag_refresh(f_msg, f_i, ref_l_os, read_l_os, LEFT);
	}
	return ret;
}*/

//Extend the frag around with f_i
//Extend PER_LEN in a while cycle, stop when get a BAD_ALIGN
//XXX If it's SHORT, extend the whole frag
//XXX Else, extend PER_LEN in a while cycle, stop when get a BAD_ALIGN
/*
int frag_extend(frag_msg *f_msg, int f_i, int FLAG, bntseq_t *bns, uint8_t *pac, char *read_seq, uint8_t *seq1, uint8_t *seq2, int seed_len)
{
	int ret, i;
	int ref_len, read_len;
	int ref_l_os, read_l_os, ref_r_os, read_r_os;
	int b_f = f_msg->fa_msg[f_i].b_f;
	int ff, m;
	char *seq_s2 = malloc(16385*sizeof(char));
	if (FLAG == RIGHT)
	{
		if (b_f == 0)
		{
			//printf("BAD\n");
			return BAD_ALIGN;
		}
		ref_len = read_len = seed_len;
		for (i = 0; i < f_msg->fa_msg[b_f-1].per_n; i++)
		{
			if (f_msg->fa_msg[f_i].srand == -1)
				pac2fa_core(bns, pac, f_msg->fa_msg[f_i].chr, f_msg->fa_msg[f_i].ex_ref_end-1-ref_len, &ref_len, f_msg->fa_msg[f_i].srand, &ff, seq1);
			else
				pac2fa_core(bns, pac, f_msg->fa_msg[f_i].chr, f_msg->fa_msg[f_i].ex_ref_end, &ref_len, f_msg->fa_msg[f_i].srand, &ff, seq1);
			strncpy(seq_s2, read_seq + (f_msg->fa_msg[b_f-1].read_begin - 1 + i * seed_len), read_len);
			for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
			ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
			frag_refresh(f_msg, f_i, ref_r_os, read_r_os, RIGHT);
			if (ret == BAD_ALIGN)
				return BAD_ALIGN;
		}
		frag_covered(f_msg, f_i, b_f-1);
	}
	else	//LEFT
	{
		if (b_f == f_msg->frag_num-1)
			return BAD_ALIGN;
		ref_len = read_len = seed_len;
		for (i = 1; i <= f_msg->fa_msg[b_f+1].per_n; i++)
		{
			if (f_msg->fa_msg[f_i].srand == -1)
				pac2fa_core(bns, pac, f_msg->fa_msg[f_i].chr, f_msg->fa_msg[f_i].ex_ref_begin, &ref_len, f_msg->fa_msg[f_i].srand, &ff, seq1);
			else
				pac2fa_core(bns, pac, f_msg->fa_msg[f_i].chr, f_msg->fa_msg[f_i].ex_ref_begin-1-ref_len, &ref_len, f_msg->fa_msg[f_i].srand, &ff, seq1);
			strncpy(seq_s2, read_seq + (f_msg->fa_msg[b_f+1].read_end - i * seed_len), read_len);
			for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
			ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);

			frag_refresh(f_msg, f_i, ref_l_os, read_l_os, LEFT);
			if (ret == BAD_ALIGN)
				return BAD_ALIGN;
		}
		frag_covered(f_msg, f_i, b_f+1);
	}
	return GOOD_ALIGN;
}
*/

//XXX combine get_ref_intv & get_read_intv
int get_ref_intv(uint8_t *seq1, int *N_flag, 
				 bntseq_t *bns, uint8_t *pac, 
				 aln_msg *a_msg, 
				 int seed1_i, int seed1_aln_i, int seed2_i, int seed2_aln_i, 
				 int seed_len)
{
	int64_t start;
	int32_t len;
	int N_len;
	//if (a_msg[seed1_i].at[seed1_aln_i].nsrand == 1)	//'+' srand
	{
		start = a_msg[seed1_i].at[seed1_aln_i].offset + seed_len - 1 + a_msg[seed1_i].at[seed1_aln_i].len_dif;
		len = a_msg[seed2_i].at[seed2_aln_i].offset - 1 - start;
		if (len <= 0) return 0;
		//printf("offset: %lld ", (long long)start);
		pac2fa_core(bns, pac, a_msg[seed1_i].at[seed1_aln_i].chr, start, &len, 1, N_flag, &N_len, seq1);
	}
	/*else	//'-' srand
	{
		start = a_msg[seed2_i].at[seed2_aln_i].offset + seed_len - 1 + a_msg[seed2_i].at[seed2_aln_i].len_dif;
		len = a_msg[seed1_i].at[seed1_aln_i].offset - 1 - start;
		if (len <= 0) return 0;
		pac2fa_core(bns, pac, a_msg[seed1_i].at[seed1_aln_i].chr, start, &len, -1, &flag, seq1);
	}*/
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
				  int seed_len)
{
	int32_t i;
	int j;

	//set band-width
	*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
	
	
	//convert char seq to int seq

	if (a_msg[seed1_i].at[seed1_aln_i].nsrand == 1)
	{
		*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
		for (j=0, i = a_msg[seed1_i].read_id * 2 * seed_len - seed_len; i < (a_msg[seed2_i].read_id-1) * 2 * seed_len; ++j, ++i)
			seq2[j] = nst_nt4_table[(int)read_seq[i]];
	}
	else
	{
		//printf("read: %d, %d ", a_msg[seed1_i].read_id,a_msg[seed1_i].read_id * 2 * seed_len);
		*band_width = 2 * (((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) << 1) - 1) * MAXOFTWO(a_msg[seed1_i].at[seed1_aln_i].bmax, a_msg[seed2_i].at[seed2_aln_i].bmax);
		for (j=0, i = a_msg[seed1_i].read_id * 2 * seed_len; i < (a_msg[seed2_i].read_id * 2 - 1) * seed_len; ++j, ++i)
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

/*int get_split_ref(uint8_t *ref1_seq, int *len1, int64_t *offset_1, uint8_t *ref2_seq, int *len2, int64_t *offset_2, bntseq_t *bns, uint8_t *pac, aln_msg *a_msg, int seed1_i, int seed1_aln_i, int seed2_i, int seed2_aln_i, int seed_len)
{
	int64_t start;
	int flag1, flag2;

	if (a_msg[seed1_i].at[seed1_aln_i].nsrand == 1)	//'+' srand
	{
		start = a_msg[seed1_i].at[seed1_aln_i].offset + seed_len - 1 + a_msg[seed1_i].at[seed1_aln_i].len_dif;
		*offset_1 = start;
		*len2 = *len1 = ((a_msg[seed2_i].read_id - a_msg[seed1_i].read_id) * 2 - 1) * seed_len;
		pac2fa_core(bns, pac, a_msg[seed1_i].at[seed1_aln_i].chr, start, len1, 1, &flag1, ref1_seq);

		start = a_msg[seed2_i].at[seed2_aln_i].offset - *len2 - 1;
		*offset_2 = start;
		pac2fa_core(bns, pac, a_msg[seed2_i].at[seed2_aln_i].chr, start, len2, 1, &flag2, ref2_seq);
	}
	else	//'-' srand
	{
		printf("XXX - srand.\n");
	}
	return 0;
}*/

void printcigar(uint32_t *cigar, int cigar_len)
{
	int i;

    //fprintf(stdout, "cigar %d: ", cigar_len);
	for (i = 0; i < cigar_len; i++)
		fprintf(stdout, "%d%c", (int)(cigar[i]>>4), "MIDNSHP=X"[(int)(cigar[i] & 0xf)]);
	//fprintf(stdout, "\t");
}

void printnst(char *msg1, uint8_t *seq, int len, char *msg2)
{
	int i;
	printf("%s", msg1);
	for (i = 0; i < len; ++i) printf("%c", "ACGTN"[(int)seq[i]]);
	printf("%s", msg2);
}
//merge adjacent same operations XXX
void _push_cigar(uint32_t **cigar, int *cigar_len, int *cigar_m, uint32_t *_cigar, int _cigar_len)
{
    if (_cigar_len == 0) return;
	int i,j;
	i = *cigar_len;
	if (((i-1) >= 0) && (((*cigar)[i-1] & 0xf) == (_cigar[0] & 0xf)))
	{
		(*cigar)[i-1] += ((_cigar[0] >> 4) << 4);
		j = 1;
	}
	else j = 0;

	for (; j < _cigar_len; ++i,++j) {
		if (i == *cigar_m) {
			(*cigar_m) <<= 1;
			(*cigar) = (uint32_t*)realloc(*cigar, (*cigar_m) * sizeof (uint32_t));
			if ((*cigar) == NULL)	{fprintf(stderr, "[frag_check] Memory is not enougy.\n");exit(-1);}
		}	
		(*cigar)[i] = _cigar[j];
	}
	*cigar_len = i;
}

//merge interval and seed to frag
//Before, frag has the first seed's msg already
void merge_cigar(frag_msg *f_msg, int f_i, 
				 uint32_t *cigar, int cigar_len, 
				 int reflen, int readlen, 
				 bntseq_t *bns, uint8_t *pac, 
				 char *read_seq)
{
	int *fc_len = &(f_msg->fa_msg[f_i].cigar_len);
	uint32_t **fcigar = &(f_msg->fa_msg[f_i].cigar);
	//printcigar(*fcigar, *fc_len);
	//printcigar(cigar, cigar_len);
	//refresh cigar msg
	//XXX
	if (((*fcigar)[*fc_len-1] & 0xf) != CMATCH || (cigar[0] & 0xf) != CMATCH)	//bound-repair
	{
		    /* seq1, ref */ /*    seq2, read     */
		int len1, len11=0,  len2, len21=0, len22=0, len_dif=0;
		int bd_cigar_len=0, b=0, min_b, score, ci=0;
		int N_flag, N_len;
		int i, j, left=1, right=1;
		uint8_t *seq1, *seq2;
		uint32_t *bd_cigar = 0;
		int64_t ref_start;
		int32_t read_start;
		//XXX 
		int md = 3;

		while (1) {	//get a right boundary-cigar(without 'I'/'D' on the head or tail)
			//printf("a");printcigar(*fcigar, *fc_len);
			//printf("b");printcigar(cigar+ci, cigar_len-ci);
			//extend 3-match boundary
			if (left) {
				while ((*fc_len) >= 1) {	// get a nM, where n > 3
					if (((*fcigar)[*fc_len-1] & 0xf) == CMATCH && ((*fcigar)[*fc_len-1]>>4) > md) {
						(*fcigar)[*fc_len-1] -= (md << 4); len21 += md; break;}
					else if (((*fcigar)[*fc_len-1] & 0xf) == CMATCH) {
						len21 += (int)((*fcigar)[*fc_len-1] >> 4); --(*fc_len); }
					else if (((*fcigar)[*fc_len-1] & 0xf) == CINS) {
						len21 += (int)((*fcigar)[*fc_len-1] >> 4); len_dif -= (int)((*fcigar)[*fc_len-1] >> 4); b += (int)((*fcigar)[*fc_len-1]>>4); --(*fc_len); }
					else if (((*fcigar)[*fc_len-1] & 0xf) == CDEL) {
						len_dif += (int)((*fcigar)[*fc_len-1] >> 4); b += (int)((*fcigar)[*fc_len-1] >> 4); --(*fc_len); }
					else { fprintf(stderr, "[merge cigar] Unexpected: "); printcigar((*fcigar), (*fc_len)); exit(-1);}
				}
				/*
				if (((*fcigar)[*fc_len-1] & 0xf) == CMATCH) {
					if ((int)((*fcigar)[*fc_len-1] >> 4) > 3) { 
						(*fcigar)[*fc_len-1] -= (0x3 << 4); len21 += 3;
					} else { 
						len21 += (int)((*fcigar)[*fc_len-1] >> 4); (*fc_len)--; }
				} else {	//CINS or CDEL
					while (1) {	//skip all the CIN-DEL opreations
						if (((*fcigar)[*fc_len-1] & 0xf) == CINS) { 
							len21 += (int)((*fcigar)[*fc_len-1]>>4); len_dif -= (int)((*fcigar)[*fc_len-1]>>4); b += (int)((*fcigar)[*fc_len-1]>>4); 
						} else { 
							len_dif += (int)((*fcigar)[*fc_len-1] >> 4); b += (int)((*fcigar)[*fc_len-1] >> 4);}
						(*fc_len)--;
						if (((*fcigar)[*fc_len-1] & 0xf) == CMATCH) break; 
						if (*fc_len == 1) {fprintf(stderr, "cin-del: "); printcigar(*fcigar, *fc_len);}
					}
					if ((int)((*fcigar)[*fc_len-1] >> 4) > 3) { 
						(*fcigar)[*fc_len-1] -= (0x3 << 4); len21 += 3;
					} else { 
						len21 += (int)((*fcigar)[*fc_len-1] >> 4); (*fc_len)--; }
				}*/
				len11 = len21 + len_dif;
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
						len22 += (int)(cigar[ci] >> 4); len_dif -= (int)(cigar[ci]>>4); b += (int)(cigar[ci]>>4); ++ci; }
					else if ((cigar[ci] & 0xf) == CDEL) {
						len_dif += (int)(cigar[ci] >> 4); b += (int)(cigar[ci]>>4); ++ci; }
					else { fprintf(stderr, "[merge cigar] Unexpected: "); printcigar((cigar+ci), cigar_len-ci); exit(-1);}
				}
				/*if ((cigar[ci] & 0xf) == CMATCH) {
					if ((int)(cigar[ci] >> 4) > 3) { 
						cigar[ci] -= (0x3 << 4); len22 += 3;
					} else { 
						len22 += (int)(cigar[ci] >> 4); ci++; }
				} else {	//CINS or CDEL
					while (1) {	//skip all the CIN-DEL opreations
						if ((cigar[ci] & 0xf) == CINS) {
							len22 += (int)(cigar[ci]>>4); len_dif -= (int)(cigar[ci]>>4); b += (int)(cigar[ci]>>4);
						} else { 	//CDEL
							len_dif += (int)(cigar[ci]>>4); b += (int)(cigar[ci]>>4);	}
						ci++;
						if ((cigar[ci] & 0xf) == CMATCH) break;
						if (ci == cigar_len-1) {fprintf(stderr, "cin-del: "); printcigar(cigar, cigar_len);}
					}
					if ((int)(cigar[ci] >> 4) > 3) {
						cigar[ci] -= (0x3 << 4); len22 += 3;
					} else {
						len22 += (int)(cigar[ci] >> 4); ci++;}
				}*/
			}
			len2 = len21+len22; len1 = len2 + len_dif;
			//printf("len1: %d len21: %d len22: %d len_dif %d len2: %d\t", len1, len21, len22, len_dif, len2);
			min_b = abs(len_dif) + md; b = MAXOFTWO(b, min_b);
			seq1 = (uint8_t *)malloc((len1+1) * sizeof(uint8_t)); seq2 = (uint8_t *)malloc((len2+1) * sizeof(uint8_t));
			pac2fa_core(bns, pac, f_msg->fa_msg[f_i].chr, ref_start-1, &len1, 1, &N_flag, &N_len, seq1);
			for (i=0, j=read_start-1; j<read_start+len2-1; ++i, ++j) seq2[i] = nst_nt4_table[(int)read_seq[j]];
			//printnst("seq1: ", seq1, len1, "\t"); printnst("seq2: ", seq2, len2, "\t");
			//fprintf(stdout, "len1: %d, len2: %d\n", len1, len2);
			//score = ksw_global(len2, seq2, len1, seq1, 5, &sc_mat, 2, 1, b, &bd_cigar_len, &bd_cigar);
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
			//else printf("go on\n");
			left = right = 1;
			//else {fprintf(stdout,"bd: %d %ld %d %d ", read_start-1, ref_start-1, len2, len1); /*printcigar(bd_cigar, bd_cigar_len);*/}
			free(bd_cigar);
		}
		//printcigar(fcigar, *fc_len);
		//printcigar(bd_cigar, bd_cigar_len);
		//printcigar(cigar+ci, cigar_len-ci);
		_push_cigar(fcigar, fc_len, &(f_msg->fa_msg[f_i].cigar_max), bd_cigar, bd_cigar_len);
		_push_cigar(fcigar, fc_len, &(f_msg->fa_msg[f_i].cigar_max), cigar+ci, cigar_len-ci);
		free(bd_cigar);
	}
	else _push_cigar(fcigar, fc_len, &(f_msg->fa_msg[f_i].cigar_max), cigar, cigar_len);
	//refresh coordinates msg
	f_msg->fa_msg[f_i].cigar_ref_end += reflen;
	f_msg->fa_msg[f_i].cigar_read_end += readlen;
	f_msg->fa_msg[f_i].len_dif += (reflen-readlen);
}

//Extend the intervals between seeds in the same frag
//XXX un-seed region: end of + srand and start of - end 
int frag_extend(frag_msg *f_msg, aln_msg *a_msg, int f_i, 
				bntseq_t *bns, uint8_t *pac, char *read_seq, 
				uint8_t *seq1, uint8_t *seq2, 
				int seed_len, 
				int **line_tri)
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
	}
	else	//'-' srand
	{
		last_i = f_msg->fa_msg[f_i].seed_i[0];
		last_aln_i = f_msg->fa_msg[f_i].seed_aln_i[0];
	}

	//copy the first seed's cigar to frag
	{
		_push_cigar(&(f_msg->fa_msg[f_i].cigar), &(f_msg->fa_msg[f_i].cigar_len), &(f_msg->fa_msg[f_i].cigar_max), a_msg[last_i].at[last_aln_i].cigar, a_msg[last_i].at[last_aln_i].cigar_len);
		/*f_msg->fa_msg[f_i].cigar_len = a_msg[last_i].at[last_aln_i].cigar_len;
		for (j = 0; j < f_msg->fa_msg[f_i].cigar_len; j++)
			f_msg->fa_msg[f_i].cigar[j] = a_msg[last_i].at[last_aln_i].cigar[j];*/
		f_msg->fa_msg[f_i].len_dif = a_msg[last_i].at[last_aln_i].len_dif;
		//start/end: 1-base
		f_msg->fa_msg[f_i].cigar_ref_start = a_msg[last_i].at[last_aln_i].offset;
		f_msg->fa_msg[f_i].cigar_ref_end = a_msg[last_i].at[last_aln_i].offset + seed_len -1 + a_msg[last_i].at[last_aln_i].len_dif;

		//read_start: '-' srand XXX
		f_msg->fa_msg[f_i].cigar_read_start = (a_msg[last_i].read_id - 1) * 2 * seed_len+1 + seed_len;		//1-based
		f_msg->fa_msg[f_i].cigar_read_end = f_msg->fa_msg[f_i].cigar_read_start - 1 + seed_len;	//1-based
		f_msg->fa_msg[f_i].edit_dis = a_msg[last_i].at[last_aln_i].edit_dis;	

		//XXX f_msg->fa_msg[f_i].cigar_ref_start = a_msg[last_i].at[last_aln_i].offset + seed_len -1 + a_msg[last_i].at[last_aln_i].len_dif;
		//_msg->fa_msg[f_i].cigar_ref_end = a_msg[last_i].at[last_aln_i].offset;
	}
	if (f_msg->fa_msg[f_i].srand == 1)
	{
		for (--i; i >= 0; --i)
		{
			seed_i = f_msg->fa_msg[f_i].seed_i[i];
			seed_aln_i = f_msg->fa_msg[f_i].seed_aln_i[i];
			//get ref seq and seed-interval seq
			if ((len1 = get_ref_intv(seq1, &N_flag, bns, pac, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, seed_len)) == 0) { fprintf(stderr, "[frag extend] frag path error : (%d,%d) -> %d %lld %d %lld\n", last_i, last_aln_i, a_msg[last_i].read_id, (long long)a_msg[last_i].at[last_aln_i].offset, a_msg[seed_i].read_id, (long long)a_msg[seed_i].at[seed_aln_i].offset); return -1;}
			len2 = get_read_intv(seq2, read_seq, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, &b, seed_len);
			min_b = abs(len2-len1)+3;
			b = MAXOFTWO(b, min_b);
			cigar_len = 0;
			uint32_t *cigar=0;
			//XXX score check
			score = ksw_global(len2, seq2, len1, seq1, 5, bwasw_sc_mat, 5, 2, b, &cigar_len, &cigar);
			merge_cigar(f_msg, f_i, cigar, cigar_len, len1, len2, bns, pac, read_seq);
			merge_cigar(f_msg, f_i, a_msg[seed_i].at[seed_aln_i].cigar, a_msg[seed_i].at[seed_aln_i].cigar_len, seed_len+a_msg[seed_i].at[seed_aln_i].len_dif, seed_len, bns, pac, read_seq);//merge seed to frag
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
			if ((len1 = get_ref_intv(seq1, &N_flag, bns, pac, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, seed_len)) == 0) { fprintf(stderr, "[frag extend] frag path error : (%d,%d) -> %d %lld %d %lld\n", last_i, last_aln_i, a_msg[last_i].read_id, (long long)a_msg[last_i].at[last_aln_i].offset, a_msg[seed_i].read_id, (long long)a_msg[seed_i].at[seed_aln_i].offset); return -1;}
			len2 = get_read_intv(seq2, read_seq, a_msg, last_i, last_aln_i, seed_i, seed_aln_i, &b, seed_len);
			min_b = abs(len2-len1)+3;
			b = MAXOFTWO(b, min_b);
			cigar_len = 0;
			uint32_t *cigar=0;
			//XXX score check
			score = ksw_global(len2, seq2, len1, seq1, 5, bwasw_sc_mat, 5, 2, b, &cigar_len, &cigar);
			//merge_cigar
			merge_cigar(f_msg, f_i, cigar, cigar_len, len1, len2, bns, pac, read_seq);
			merge_cigar(f_msg, f_i, a_msg[seed_i].at[seed_aln_i].cigar, a_msg[seed_i].at[seed_aln_i].cigar_len, seed_len+a_msg[seed_i].at[seed_aln_i].len_dif, seed_len, bns, pac, read_seq);//merge seed to frag

			last_i = seed_i;
			last_aln_i = seed_aln_i;
			free(cigar);
		}
	}
	
	return 0;
}

/*void split_mapping(bntseq_t *bns, uint8_t *pac, char *read_seq, frag_msg *f_msg, aln_msg *a_msg, int32_t **hash_num, int32_t ***hash_node, hash_aln_msg **ha_msg, int *ha_num, int f1_i, int f2_i, int seed_len)
{
    int i;
	int s1_i, s1_aln_i, seed_num, s2_i, s2_aln_i;

	int split_read_len, split_ref1_len, split_ref2_len;
	uint8_t *split_read_seq, *split_ref1_seq, *split_ref2_seq;
	int64_t offset_1, offset_2;

	s1_i = f_msg->fa_msg[f1_i].seed_i[0];
	s1_aln_i = f_msg->fa_msg[f1_i].seed_aln_i[0];

	seed_num = f_msg->fa_msg[f2_i].seed_num;
	s2_i = f_msg->fa_msg[f2_i].seed_i[seed_num-1];
	s2_aln_i = f_msg->fa_msg[f2_i].seed_aln_i[seed_num-1];

	split_read_len = (2*(a_msg[s2_i].read_id-a_msg[s1_i].read_id)-1) * seed_len;
	split_ref1_len = split_read_len;
	split_ref2_len = split_read_len;

	split_read_seq = (uint8_t*)malloc(split_read_len * sizeof(uint8_t));
	split_ref1_seq = (uint8_t*)malloc(split_ref1_len * sizeof(uint8_t));
	split_ref2_seq = (uint8_t*)malloc(split_ref2_len * sizeof(uint8_t));

	get_split_ref(split_ref1_seq, &split_ref1_len, &offset_1, split_ref2_seq, &split_ref2_len, &offset_2, bns, pac, a_msg, s1_i, s1_aln_i, s2_i, s2_aln_i, seed_len);
    for (i = 0; i < split_read_len; ++i)
        split_read_seq[i] = nst_nt4_table[(int)read_seq[(a_msg[s1_i].read_id * 2 - 1)*seed_len+i]];
	//strncpy(split_read_seq, read_seq+(a_msg[s1_i].read_id * 2 - 1)*seed_len, split_read_len);

	//check SV-type
	if (a_msg[s1_i].at[s1_aln_i].chr == a_msg[s2_i].at[s2_aln_i].chr &&
		a_msg[s1_i].at[s1_aln_i].nsrand == a_msg[s2_i].at[s2_aln_i].nsrand)
	{
		int pos = (a_msg[s1_i].at[s1_aln_i].nsrand == 1?a_msg[s1_i].at[s1_aln_i].offset+seed_len-1+a_msg[s1_i].at[s1_aln_i].len_dif:
				                                        a_msg[s2_i].at[s2_aln_i].offset+seed_len-1+a_msg[s2_i].at[s2_aln_i].len_dif);
		int exp = a_msg[s1_i].at[s1_aln_i].offset + a_msg[s1_i].at[s1_aln_i].nsrand * (a_msg[s2_i].read_id - a_msg[s1_i].read_id) * 2 * seed_len;	
		int act = a_msg[s2_i].at[s2_aln_i].offset;
		int dis = a_msg[s1_i].at[s1_aln_i].nsrand * (act-exp) - ((a_msg[s1_i].at[s1_aln_i].nsrand==1)?a_msg[s1_i].at[s1_aln_i].len_dif:a_msg[s2_i].at[s2_aln_i].len_dif);

		if (dis > 0) fprintf(stdout, "%d\t%d\t%d\tDEL\n", a_msg[s1_i].at[s1_aln_i].chr, pos, dis);
		else if (dis < 0) fprintf(stdout, "%d\t%d\t%d\tINS\n", a_msg[s1_i].at[s1_aln_i].chr, pos, (0-dis));
		else fprintf(stderr, "[split-map] frag error: %d %d %lld.\n", a_msg[s1_i].at[s1_aln_i].chr, a_msg[s1_i].at[s1_aln_i].nsrand, (long long)a_msg[s1_i].at[s1_aln_i].offset);
		//split hash-map
		split_map(split_read_seq, split_read_len, split_ref1_seq, split_ref1_len, offset_1, split_ref2_seq, split_ref2_len, offset_2, 10, hash_num, hash_node, ha_msg, ha_num, 2, 16);

	}
	else	//dif srand
	{
        int pos1, pos2;
        if (a_msg[s1_i].at[s1_aln_i].nsrand == 1)
            pos1 = a_msg[s1_i].at[s1_aln_i].offset+seed_len-1;
        else pos1 = a_msg[s1_i].at[s1_aln_i].offset;
        if (a_msg[s2_i].at[s2_aln_i].nsrand == 1)
            pos2 = a_msg[s2_i].at[s2_aln_i].offset;
        else pos2 = a_msg[s2_i].at[s2_aln_i].offset+seed_len-1;
		fprintf(stdout, "%d\t%d\tcomplex SV\n", a_msg[s1_i].at[s1_aln_i].chr, pos1);
		fprintf(stdout, "%d\t%d\tcomplex SV\n", a_msg[s2_i].at[s2_aln_i].chr, pos2);
	}
	free(split_read_seq); free(split_ref1_seq); free(split_ref2_seq);
}*/

//for '-' srand:
//1. convert read seq into reverse-complement seq.
//2. split-map the re-co seq to the ref.
//3. combine with aln results of other regions.
//@func: deal with seqs between 2 frags
void split_mapping(uint32_t **split_cigar, int *split_len, int *split_m,
				   bntseq_t *bns, uint8_t *pac, 
				   char *read_seq, 
				   frag_msg *f_msg, aln_msg *a_msg, 
				   uint32_t **hash_num, uint64_t ***hash_node, 
				   int f1_i, int f2_i, 
				   int seed_len, 
				   int **line_tri)
{
    int i;
    int nsrand, s1_i, s1_aln_i, seed_num, s2_i, s2_aln_i;

	nsrand = f_msg->fa_msg[f1_i].srand;
	if (nsrand == 1)
	{
		s1_i = f_msg->fa_msg[f1_i].seed_i[0];
		s1_aln_i = f_msg->fa_msg[f1_i].seed_aln_i[0];

		seed_num = f_msg->fa_msg[f2_i].seed_num;
		s2_i = f_msg->fa_msg[f2_i].seed_i[seed_num-1];
		s2_aln_i = f_msg->fa_msg[f2_i].seed_aln_i[seed_num-1];
	}
	else
	{
		seed_num = f_msg->fa_msg[f1_i].seed_num;
		s1_i = f_msg->fa_msg[f1_i].seed_i[seed_num-1];
		s1_aln_i = f_msg->fa_msg[f1_i].seed_aln_i[seed_num-1];

		s2_i = f_msg->fa_msg[f2_i].seed_i[0];
		s2_aln_i = f_msg->fa_msg[f2_i].seed_aln_i[0];
	}

	int split_read_len, split_ref_len;
	uint8_t *split_read_seq, *split_ref_seq;
	int64_t ref_offset;
	//init value for hash-map
	int hash_len = HASH_LEN;
	int nt_n = NT_N;
	int key_len = 2;
	int hash_size = (int)pow(nt_n, key_len);
	int hash_step = HASH_STEP;
	int N_flag, N_len;

	split_read_len = (2*(a_msg[s2_i].read_id-a_msg[s1_i].read_id)-1) * seed_len;
	if (split_read_len > (int)pow((double)(2.0),(double)(16.0))) fprintf(stderr, "[split map] ERROR: split read is too long: %d\n", split_read_len);
	split_read_seq = (uint8_t*)malloc(split_read_len * sizeof(uint8_t));

	if (nsrand == 1) for (i = 0; i < split_read_len; ++i) split_read_seq[i] = nst_nt4_table[(int)read_seq[(a_msg[s1_i].read_id * 2 - 1)*seed_len+i]];
    else for (i = 0; i < split_read_len; ++i) split_read_seq[i] = nst_nt4_table[(int)read_seq[(a_msg[s1_i].read_id * 2) * seed_len + i]];//'-' srand	

	//check SV-type
	if (a_msg[s1_i].at[s1_aln_i].chr == a_msg[s2_i].at[s2_aln_i].chr &&
		a_msg[s1_i].at[s1_aln_i].nsrand == a_msg[s2_i].at[s2_aln_i].nsrand)
	{
		int64_t pos = a_msg[s1_i].at[s1_aln_i].offset+seed_len-1+a_msg[s1_i].at[s1_aln_i].len_dif;
		int64_t exp = a_msg[s1_i].at[s1_aln_i].offset + (a_msg[s2_i].read_id - a_msg[s1_i].read_id) * 2 * seed_len;	
		int64_t act = a_msg[s2_i].at[s2_aln_i].offset;
		int dis = (act-exp) - a_msg[s1_i].at[s1_aln_i].len_dif;

		//fprintf(stdout, "s1: %lld, exp_s2: %lld act_s2: %lld %d\n", (long long)pos, (long long)exp, (long long)act, dis);
		if (dis > 10)
		{	
#ifdef __DEBUG__
			fprintf(stdout, "%d\t%lld\t%d\tDEL\n", a_msg[s1_i].at[s1_aln_i].chr, (long long)pos, dis);
#endif
			//long-del XXX
			//if (dis < split_read_len)
			{
				split_ref_len = split_read_len + dis;
				split_ref_seq = (uint8_t*)malloc(split_ref_len * sizeof(uint8_t));

				//ref_offset : left most postion of split-ref(1-base)
				ref_offset = a_msg[s1_i].at[s1_aln_i].offset + seed_len + a_msg[s1_i].at[s1_aln_i].len_dif;
				//XXX s1/s2???
				//else ref_offset = a_msg[s2_i].at[s2_aln_i].offset + seed_len + a_msg[s2_i].at[s2_aln_i].len_dif; 
				pac2fa_core(bns, pac, a_msg[s1_i].at[s1_aln_i].chr, ref_offset-1, &split_ref_len, 1/*+-*/, &N_flag, &N_len, split_ref_seq);
				split_delete_map(split_cigar, split_len, split_m, split_read_seq, split_read_len, split_ref_seq, split_ref_len, ref_offset, hash_len, hash_step, hash_num, hash_node, key_len, hash_size);
			}
			//else	//dis > split_read_len
		}
		else if (dis < -10) 
		{
#ifdef __DEBUG__
			fprintf(stdout, "%d\t%lld\t%d\tINS\n", a_msg[s1_i].at[s1_aln_i].chr, (long long)pos, (0-dis));
#endif
			split_ref_len = split_read_len + dis;
			if (split_ref_len < hash_len) 
			{
				hash_len = split_ref_len;	//XXX
				if (hash_len <= key_len)
				{
					fprintf(stderr, "[frag check] split ref is too short(<= key_len, 2).\n");
					exit(-1);
				}
			}
			if (split_ref_len <= 0) {fprintf(stderr, "[split map] split length error.\n"); exit(-1);}
			split_ref_seq = (uint8_t*)malloc(split_ref_len * sizeof(uint8_t));
			//if (nsrand == 1) 
			ref_offset = a_msg[s1_i].at[s1_aln_i].offset + seed_len + a_msg[s1_i].at[s1_aln_i].len_dif;
			//else ref_offset = a_msg[s2_i].at[s2_aln_i].offset + seed_len + a_msg[s2_i].at[s2_aln_i].len_dif;
			pac2fa_core(bns, pac, a_msg[s1_i].at[s1_aln_i].chr, ref_offset-1, &split_ref_len, 1, &N_flag, &N_len, split_ref_seq);
			split_insert_map(split_cigar, split_len, split_m, split_read_seq, split_read_len, split_ref_seq, split_ref_len, ref_offset, hash_len, hash_step, hash_num, hash_node, key_len, hash_size);
		}
		//for MIS-MATCH
		else 
		{
			//XXX	MIS-MATCH
			split_ref_len = split_read_len + dis;
			split_ref_seq = (uint8_t*)malloc(split_ref_len * sizeof(uint8_t));
			ref_offset = a_msg[s1_i].at[s1_aln_i].offset + seed_len + a_msg[s1_i].at[s1_aln_i].len_dif;
			pac2fa_core(bns, pac, a_msg[s1_i].at[s1_aln_i].chr, ref_offset-1, &split_ref_len, 1, &N_flag, &N_len, split_ref_seq);
			int min_b = abs(dis)+3;
			uint32_t *cigar=0;
			int score = ksw_global(split_read_len, split_read_seq, split_ref_len, split_ref_seq, 5, bwasw_sc_mat, 5, 2, min_b, split_len, &cigar);
			//XXX 'spurious MATCH'
			if (score < 0) 
			{
				//fprintf(stdout, "cigar:\t"); printcigar(cigar, *split_len); fprintf(stdout, "band: %d\n", min_b);
				//fprintf(stdout, "ref:\t"); for (i = 0; i < split_ref_len; ++i) fprintf(stdout, "%c", "ACGT"[split_ref_seq[i]]); fprintf(stdout, "\n");
				//fprintf(stdout, "read:\t");for (i = 0; i < split_read_len; ++i) fprintf(stdout, "%c", "ACGT"[split_read_seq[i]]); *split_len = 0; fprintf(stdout, "\n");

				if (hash_both_bound_map(split_cigar, split_len, split_m, split_ref_seq, split_ref_len, split_read_seq, split_read_len, hash_num, hash_node, HASH_LEN, 2, HASH_STEP) == 1)//XXX pull trigger
				{
					if (nsrand == 1)
					{
						for(i = 0; i < f_msg->fa_msg[f2_i].pre_trigger_n; ++i)
							(*line_tri)[f_msg->fa_msg[f2_i].trigger[f_msg->fa_msg[f2_i].next_trigger_n+i]] = 1;
					}
					else
					{
						for(i = 0; i < f_msg->fa_msg[f1_i].pre_trigger_n; ++i)
							(*line_tri)[f_msg->fa_msg[f2_i].trigger[f_msg->fa_msg[f2_i].next_trigger_n+i]] = 1;
					}
				}
				//fprintf(stdout, "both: "); printcigar(*split_cigar, *split_len); fprintf(stdout, "\n");

				//cut cigar, fix right bound of 1'st half read
				/*hash_right_bound_map(&b_cigar, &cigar_len, split_ref_seq, split_ref_len, split_read_seq, split_read_len, hash_num, hash_node, HASH_LEN, 2, HASH_STEP);
				for (i = 0; i < cigar_len; ++i)
				{
					if ((b_cigar[i] & 0xf) != CSOFT_CLIP)
					{
						(*split_cigar)[*split_len] = b_cigar[i];
						(*split_len)++;
					}
				}
				fprintf(stdout, "1st: "); printcigar(b_cigar, cigar_len); fprintf(stdout, "\n");
				(*split_cigar)[*split_len] = CSOFT_CLIP;
				(*split_len)++;

				//fix left bound of 2'nd half read
				hash_left_bound_map(&b_cigar, &cigar_len, split_ref_seq, split_ref_len, split_read_seq, split_read_len, hash_num, hash_node, HASH_LEN, 2, HASH_STEP);
				for (i = 0; i < cigar_len; ++i)
				{
					if ((b_cigar[i] & 0xf) != CSOFT_CLIP)
					{
						(*split_cigar)[*split_len] = b_cigar[i];
						(*split_len)++;
					}
				}
				fprintf(stdout, "2nd: "); printcigar(b_cigar, cigar_len); fprintf(stdout, "\n");*/
			}
			else
			{
				for (i = 0; i < *split_len; ++i)
					(*split_cigar)[i] = cigar[i];
			}
			free(cigar);
		}
	}
	else	//dif srand
	{
        int pos1, pos2;
        if (a_msg[s1_i].at[s1_aln_i].nsrand == 1)
            pos1 = a_msg[s1_i].at[s1_aln_i].offset+seed_len-1;
        else pos1 = a_msg[s1_i].at[s1_aln_i].offset;
        if (a_msg[s2_i].at[s2_aln_i].nsrand == 1)
            pos2 = a_msg[s2_i].at[s2_aln_i].offset;
        else pos2 = a_msg[s2_i].at[s2_aln_i].offset+seed_len-1;
		fprintf(stderr, "%d\t%d\tcomplex SV\n", a_msg[s1_i].at[s1_aln_i].chr, pos1);
		fprintf(stderr, "%d\t%d\tcomplex SV\n", a_msg[s2_i].at[s2_aln_i].chr, pos2);
	}
	free(split_read_seq); free(split_ref_seq);
}

int check_cigar(uint32_t *cigar, int cigar_len)
{
	int i;
	int mid[5] = {0, 0, 0, 0, 0};
	for (i = 0; i < cigar_len; ++i)
		mid[(int)(cigar[i]&0xf)] += (cigar[i]>>4);
	fprintf(stdout, "\t%d M, %d I, %d D, %d S\t", mid[0], mid[1], mid[2], mid[4]);
	return 0;
}

//before first frag
//do hash map	XXX
int frag_head_bound_fix(frag_msg *f_msg, aln_msg *a_msg, 
						uint32_t **cigar, int *cigar_len, int *cigar_m,
						uint64_t *offset, bntseq_t *bns, uint8_t *pac, 
						char *read_seq, uint8_t *seq1, uint8_t *seq2, 
						int seed_len, 
						uint32_t **hash_num, uint64_t ***hash_node, 
						int **line_tri)
{
	int left_bound = f_msg->frag_left_bound;
	int frag_i, seed_x, seed_i, aln_i, i;
	int N_flag, N_len;
	int read_len, ref_len, read_start/*0-base*/;
	uint64_t ref_start;
	if (f_msg->fa_msg[0].srand == 1)	//'+' srand
	{
		frag_i = f_msg->frag_num-1;
		seed_x = f_msg->fa_msg[frag_i].seed_num - 1;
		seed_i = f_msg->fa_msg[frag_i].seed_i[seed_x];
		aln_i = f_msg->fa_msg[frag_i].seed_aln_i[seed_x];
		if (a_msg[seed_i].read_id != 1)
		//if ((a_msg[seed_i].read_id - left_bound) > 1)
		{
			read_len = ((left_bound==0)?0:seed_len) + (a_msg[seed_i].read_id - left_bound - 1) * 2 * seed_len;
			read_start = (left_bound==0)?0:((left_bound * 2 - 1 ) * seed_len);
		}
		else
		{
			(*cigar_len) = 0;
			*offset = a_msg[seed_i].at[aln_i].offset;
			return 0;
		}
	}
	else	//'-' srand
	{
		frag_i = 0; seed_x = 0;
		seed_i = f_msg->fa_msg[0].seed_i[0];
		aln_i = f_msg->fa_msg[0].seed_aln_i[0];
		read_len = ((left_bound==0)?f_msg->last_len:seed_len) + (a_msg[f_msg->fa_msg[0].seed_i[0]].read_id - 1 - left_bound) * 2 * seed_len;
		read_start= (left_bound==0)?0:(f_msg->last_len+seed_len+(left_bound-1)*2*seed_len);
	}
	//hash map
	{	
		//XXX un-overlap
		int hash_len = HASH_LEN, hash_step = HASH_STEP, hash_key = 2;
		(*offset) = a_msg[seed_i].at[aln_i].offset;
		ref_len = read_len + hash_step * 2;	//XXX
		//read	XXX check for 'N'
		for (i = 0; i < read_len; ++i)
			seq1[i] = nst_nt4_table[(int)read_seq[read_start+i]];
		//ref
		ref_start = a_msg[seed_i].at[aln_i].offset - ref_len;	//1-base
		//printf("head, ref-start: %lld, read-start: %d\n", ref_start, read_start);
		//ref XXX check for 'N'
		pac2fa_core(bns, pac, a_msg[seed_i].at[aln_i].chr, ref_start-1/*0-base*/, &ref_len, 1, &N_flag, &N_len, seq2);


		//printf("read: %s\nref: %s\n", re_read, re_ref);
		if (hash_left_bound_map(cigar, cigar_len, cigar_m, seq2, ref_len, seq1, read_len, hash_num, hash_node, hash_len, hash_key, hash_step) == 1)
		{
			if (f_msg->fa_msg[0].srand == 1)	//'+' srand
			{
				for (i = 0; i < f_msg->fa_msg[frag_i].pre_trigger_n; ++i)
					(*line_tri)[f_msg->fa_msg[frag_i].trigger[f_msg->fa_msg[frag_i].next_trigger_n+i]] = 1;
			}
			else
			{
				for (i = 0; i < f_msg->fa_msg[0].next_trigger_n; ++i)
					(*line_tri)[f_msg->fa_msg[0].trigger[i]] = 1;
			}
		}
		//*offset XXX
		for (i = 0; i < *cigar_len; ++i)
		{
			if (((*cigar)[i] & 0xf) == CMATCH || ((*cigar)[i] & 0xf) == CDEL)
				(*offset) -= ((*cigar)[i]>>4);
		}
	}

	return 0;
}

int frag_tail_bound_fix(frag_msg *f_msg, aln_msg *a_msg, 
						uint32_t **cigar, int *cigar_len, int *cigar_m,
						bntseq_t *bns, uint8_t *pac, 
						char *read_seq, uint8_t *seq1, uint8_t *seq2, 
						int seed_len, 
						uint32_t **hash_num, uint64_t ***hash_node, 
						int **line_tri)
{
	int right_bound = f_msg->frag_right_bound;
	int frag_i, seed_x, seed_i, aln_i, i;
	int N_flag, N_len;
	int read_len, read_start, ref_len;
	uint64_t ref_start;

	if (f_msg->fa_msg[0].srand == 1)	//'+' srand
	{
		//frag_i = 0; seed_x = 0;
		seed_i = f_msg->fa_msg[0].seed_i[0];
		aln_i = f_msg->fa_msg[0].seed_aln_i[0];
		read_start = ((a_msg[seed_i].read_id - 1) * 2 + 1) * seed_len;
		read_len = ((right_bound==f_msg->seed_all+1)?f_msg->last_len:seed_len) + (right_bound - 1 - a_msg[seed_i].read_id) * 2 * seed_len;
	}
	else
	{
		frag_i = f_msg->frag_num - 1;
		seed_x = f_msg->fa_msg[frag_i].seed_num - 1;
		seed_i = f_msg->fa_msg[frag_i].seed_i[seed_x];
		aln_i = f_msg->fa_msg[frag_i].seed_aln_i[seed_x];
		if (a_msg[seed_i].read_id != f_msg->seed_all)
		{
			read_start = a_msg[seed_i].read_id * 2 * seed_len;
			read_len = ((right_bound== f_msg->seed_all+1)?0:seed_len) + (right_bound - 1 - a_msg[seed_i].read_id) * 2 * seed_len;
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
		int hash_len = HASH_LEN, hash_step = HASH_STEP, hash_key = 2;
		//if (read_len > 20000){ (*cigar_len) = 0;  fprintf(stderr, "long tail boundary.\t"); return 0;}
		ref_len = read_len + hash_step * 2;
		//read
		for (i = 0; i < read_len; ++i)
			seq1[i] = nst_nt4_table[(int)read_seq[read_start+i]];
		//ref
		ref_start = a_msg[seed_i].at[aln_i].offset + seed_len;	//1-base
		//printf("tail, ref-start: %lld, read-start: %d\n", ref_start, read_start);
		pac2fa_core(bns, pac, a_msg[seed_i].at[aln_i].chr, ref_start-1/*0-base*/, &ref_len, 1, &N_flag, &N_len, seq2);

		//fprintf(stderr, "read: %s\nref: %s\n", read, ref);
		if (hash_right_bound_map(cigar, cigar_len, cigar_m, seq2, ref_len, seq1, read_len, hash_num, hash_node, hash_len, hash_key, hash_step) == 1)	//XXX pull trigger
		{
			if (f_msg->fa_msg[0].srand == 1)	//'+' srand
			{
				for (i = 0; i < f_msg->fa_msg[0].next_trigger_n; ++i)
					(*line_tri)[f_msg->fa_msg[0].trigger[i]] = 1;
			}
			else
			{
				for (i = 0; i < f_msg->fa_msg[frag_i].pre_trigger_n; ++i)
					(*line_tri)[f_msg->fa_msg[frag_i].trigger[f_msg->fa_msg[frag_i].next_trigger_n+i]] = 1;
			}
		}
	}

	return 0;
}

//read_seq: char or uint8_t?
int frag_check(char *read_name, bntseq_t *bns, uint8_t *pac, const char *read_prefix, 
			   char *read_seq, int read_len, 
			   frag_msg **f_msg, int line_n, int *line_tri, 
			   aln_msg *a_msg, 
			   uint32_t **hash_num, uint64_t ***hash_node, 
			   int seed_len)
{
	//if (strcmp(read_name, "donor1_17444953_17544952_1:0:0_1049:0:0_33/1_1") ==0)// && strcmp(read_name, "donor1_247131951_247231950_0:0:0_989:0:0_f1/1_1")!=0)
	//	printf("debug");
	//	return 0;
	int i, j;
	//XXX len
	int max_len = 100000;
	uint8_t *seq1 = (uint8_t*)malloc((max_len+1)*sizeof(uint8_t));
	uint8_t *seq2 = (uint8_t*)malloc((max_len+1)*sizeof(uint8_t));
	
	//initialization of aln-res
	aln_res a_res;
	a_res.res_m = line_n;
	a_res.cur_res_n = 0;
	a_res.la = (line_aln_res*)malloc(line_n * sizeof(line_aln_res));
	for (i = 0; i < line_n; ++i)
	{
		a_res.la[i].c_m = 100;
		a_res.la[i].cigar = (uint32_t*)malloc((a_res.la[0].c_m) * sizeof(uint32_t));
	}
	
	uint32_t *cigar;
	int cigar_len, cigar_m;
	uint64_t offset;
	cigar_m = 100;
	cigar = (uint32_t*)malloc(cigar_m * sizeof(uint32_t));
	read_len = 100000;
	// check for every line
	for (j = 0; j < line_n; ++j)
	{
        if (line_tri[j] == 0) continue;
		//res	XXX for multi-line
		a_res.la[a_res.cur_res_n].cigar_len = 0;
		a_res.la[a_res.cur_res_n].nsrand = (((*f_msg)+j)->fa_msg[0].srand == 1)?1:0;
		a_res.la[a_res.cur_res_n].chr = ((*f_msg)+j)->fa_msg[0].chr;
		//a_res.offset = 0;
		// extend once
		// for every frag : take it as a RIGHT frag
		// for every interval : SV breakpoint(s).
		if (((*f_msg)+j)->fa_msg[0].srand == 1)
		{
			//head 'S'
				if (((*f_msg)+j)->frag_left_bound > 0)
				{
					cigar_len = 1;
					cigar[0] = (((((*f_msg)+j)->frag_left_bound * 2 - 1) * seed_len) << 4) | CSOFT_CLIP;
					_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
				}
			//fix the boundary blank before first frag, if it is existed
			frag_head_bound_fix((*f_msg)+j, a_msg, &cigar, &cigar_len, &cigar_m, &offset, bns, pac, read_seq, seq1, seq2, seed_len, hash_num, hash_node, &line_tri);
			a_res.la[a_res.cur_res_n].offset = offset;
			//if (cigar_len != 0) {printf("head:\t"); printcigar(cigar, cigar_len); printf("\t");}
			_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
			for (i = ((*f_msg)+j)->frag_num-1; i > 0; --i)
			{
				frag_extend((*f_msg)+j, a_msg, i, bns, pac, read_seq, seq1, seq2, seed_len, &line_tri);
				_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), 
						&(a_res.la[a_res.cur_res_n].c_m), ((*f_msg)+j)->fa_msg[i].cigar, ((*f_msg)+j)->fa_msg[i].cigar_len);
				split_mapping(&cigar, &cigar_len, &cigar_m, bns, pac, read_seq, ((*f_msg)+j), a_msg, hash_num, hash_node, i, i-1, seed_len, &line_tri);
				_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
			}
			frag_extend((*f_msg)+j, a_msg, i, bns, pac, read_seq, seq1, seq2, seed_len, &line_tri);
			_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), ((*f_msg)+j)->fa_msg[i].cigar, ((*f_msg)+j)->fa_msg[i].cigar_len);
			//fix the boundary blank after the last frag
			frag_tail_bound_fix((*f_msg)+j, a_msg, &cigar, &cigar_len, &cigar_m, bns, pac, read_seq, seq1, seq2, seed_len, hash_num, hash_node, &line_tri);
			//if (cigar_len != 0) {printf("tail:\t"); printcigar(cigar, cigar_len); printf("\t");}
			_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
			//tail 'S'
				if (((*f_msg)+j)->frag_right_bound <= ((*f_msg)+j)->seed_all)
				{
					cigar_len = 1;
					cigar[0] = (((((*f_msg)+j)->seed_all - ((*f_msg)+j)->frag_right_bound + 1) * 2 * seed_len - seed_len + ((*f_msg)+j)->last_len) << 4) | CSOFT_CLIP;
					_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
				}
		}
		else	//'-' srand
		{
			//convert into rev-com
			char *reco_read_seq = (char*)malloc(read_len * sizeof(char));
			for (i = 0; i < read_len; ++i) reco_read_seq[i] = (read_seq[read_len-1-i]=='A')?'T':((read_seq[read_len-1-i]=='C')?'G':((read_seq[read_len-1-i]=='G')?'C':((read_seq[read_len-1-i]=='T')?'A':'N')));
			for (i = 0; i < ((*f_msg)+j)->seed_all; ++i) a_msg[i].read_id = (((*f_msg)+j)->seed_all + 1 - a_msg[i].read_id);
			int tmp = ((*f_msg)+j)->frag_left_bound;
			((*f_msg)+j)->frag_left_bound = ((*f_msg)+j)->seed_all + 1 - ((*f_msg)+j)->frag_right_bound;
			((*f_msg)+j)->frag_right_bound = ((*f_msg)+j)->seed_all + 1 - tmp;

			//head 'S'
				if (((*f_msg)+j)->frag_left_bound > 0)
				{
					cigar_len = 1;
					cigar[0] = ((((*f_msg)+j)->frag_left_bound * 2 * seed_len - seed_len + ((*f_msg)+j)->last_len) << 4) | CSOFT_CLIP;
					_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
				}
			//fix the boundary blank before first frag
			frag_head_bound_fix((*f_msg)+j, a_msg, &cigar, &cigar_len, &cigar_m, &offset, bns, pac, reco_read_seq, seq1, seq2, seed_len, hash_num, hash_node, &line_tri);
			//if (cigar_len != 0) {printf("head:\t"); printcigar(cigar, cigar_len); printf("\t");}
			a_res.la[a_res.cur_res_n].offset = offset;
			_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
			for (i = 0; i < ((*f_msg)+j)->frag_num-1; ++i)
			{
				frag_extend((*f_msg)+j, a_msg, i, bns, pac, reco_read_seq, seq1, seq2, seed_len, &line_tri);
				_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), ((*f_msg)+j)->fa_msg[i].cigar, ((*f_msg)+j)->fa_msg[i].cigar_len);
				split_mapping(&cigar, &cigar_len, &cigar_m, bns, pac, reco_read_seq, (*f_msg)+j, a_msg, hash_num, hash_node, i, i+1, seed_len, &line_tri);
				_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
			}
			frag_extend(((*f_msg)+j), a_msg, i, bns, pac, reco_read_seq, seq1, seq2, seed_len, &line_tri);
			_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), ((*f_msg)+j)->fa_msg[i].cigar, ((*f_msg)+j)->fa_msg[i].cigar_len);
			//fix the boundary blank after the last frag
			frag_tail_bound_fix(((*f_msg)+j), a_msg, &cigar, &cigar_len, &cigar_m, bns, pac, reco_read_seq, seq1, seq2, seed_len, hash_num, hash_node, &line_tri);
			//if (cigar_len != 0) {printf("tail:\t"); printcigar(cigar, cigar_len); printf("\t");}
			_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
			//tail 'S'
				if (((*f_msg)+j)->frag_right_bound <= ((*f_msg)+j)->seed_all)
				{
					cigar_len = 1;
					cigar[0] = ((((((*f_msg)+j)->seed_all - ((*f_msg)+j)->frag_right_bound + 1) * 2 - 1) * seed_len) << 4) | CSOFT_CLIP;
					_push_cigar(&(a_res.la[a_res.cur_res_n].cigar), &(a_res.la[a_res.cur_res_n].cigar_len), &(a_res.la[a_res.cur_res_n].c_m), cigar, cigar_len);
				}
			free(reco_read_seq);
			for (i = 0; i < ((*f_msg)+j)->seed_all; ++i) a_msg[i].read_id = (((*f_msg)+j)->seed_all + 1 - a_msg[i].read_id);
		}
		
		//fprintf(stdout, "%s\t%d\t%c\t%lld\t", read_name, a_res.chr, "-++"[a_res.nsrand], (long long)a_res.offset); printcigar(a_res.cigar, a_res.cigar_len); check_cigar(a_res.cigar, a_res.cigar_len); 
		//fprintf(stdout, "\n"); 
		++a_res.cur_res_n;
	}

	//fprintf(stdout, "main\t");
	//if (a_res.cur_res_n > 1) fprintf(stdout, "multi\t");
	for (i = 0; i < a_res.cur_res_n; ++i)
	{
		fprintf(stdout, "%s\t%d\t%c\t%lld\t", read_name, a_res.la[i].chr, "-++"[a_res.la[i].nsrand], (long long)a_res.la[i].offset); printcigar(a_res.la[i].cigar, a_res.la[i].cigar_len); check_cigar(a_res.la[i].cigar, a_res.la[i].cigar_len); 
		fprintf(stdout, "\n"); 
	}


	free(seq1); free(seq2); free(cigar);
	for (i = 0; i < a_res.res_m; ++i) free(a_res.la[i].cigar); free(a_res.la);
	return 0;
}
