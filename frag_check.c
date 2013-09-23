/*
 * frag_check.c
 *
 * function: transform raw frag.msg to detailed frag.msg
 * method:   extend-ssw
 *
 * Created by Yan Gao on 08/29/2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include "frag_check.h"
#include "bntseq.h"
#include "ssw.h"
#include "kseq.h"

#define PER_LEN 100

KSEQ_INIT(gzFile, gzread)


/*
 * parameter: frag.msg
 */

int extend_ssw(int8_t* ref_seq, int8_t* read_seq, int ref_len, int read_len, int* ref_l_os, int* read_l_os, int* ref_r_os, int* read_r_os);

frag_msg* frag_init_msg(int frag_max)
{
	frag_msg *f_msg = (frag_msg*)malloc(sizeof(frag_msg));

	f_msg->frag_max = frag_max; 
	f_msg->frag_num = 0;
	f_msg->last_len = PER_LEN;	//XXX
	f_msg->seed_num = 500;		//XXX
	f_msg->chr= (int *)malloc(frag_max * sizeof(int));
	f_msg->srand = (int *)malloc(frag_max * sizeof(int));

	f_msg->read_begin = (int *)malloc(frag_max * sizeof(int));
	f_msg->read_end = (int *)malloc(frag_max * sizeof(int));
	f_msg->ref_begin = (int *)malloc(frag_max * sizeof(int));
	f_msg->ref_end = (int *)malloc(frag_max * sizeof(int));

	f_msg->ex_read_begin = (int *)malloc(frag_max * sizeof(int));
	f_msg->ex_read_end = (int *)malloc(frag_max * sizeof(int));
	f_msg->ex_ref_begin = (int *)malloc(frag_max * sizeof(int));
	f_msg->ex_ref_end = (int *)malloc(frag_max * sizeof(int));

	f_msg->per_n = (int *)malloc(frag_max * sizeof(int));

	f_msg->flag = (int *)malloc(frag_max * sizeof(int));
	f_msg->b_f = (int *)malloc(frag_max * sizeof(int));

	return f_msg;
}

void frag_free_msg(frag_msg *f_msg)
{
	free(f_msg->chr);
	free(f_msg->srand);
	free(f_msg->ref_begin);
	free(f_msg->ref_end);
	free(f_msg->read_begin);
	free(f_msg->read_end);
	free(f_msg->ex_ref_begin);
	free(f_msg->ex_ref_end);
	free(f_msg->ex_read_end);
	free(f_msg->ex_read_begin);
	free(f_msg->per_n);
	free(f_msg->flag);
	free(f_msg->b_f);

	free(f_msg);
}

int frag_set_msg(aln_msg *a_msg, int seed_i, int aln_i, int FLAG, frag_msg *f_msg, int frag_i, int seed_len)//FLAG 0: start/1:end
{
	if (FLAG == 1)	//end
	{
		f_msg->chr[frag_i] = a_msg[seed_i].chr[aln_i];		//chr#
		f_msg->srand[frag_i] = a_msg[seed_i].nsrand[aln_i];	//+:1/-:-1
		//f_msg->ref_2[frag_i] = a_msg[seed_i].offset[aln_i];
		//f_msg->read_2[frag_i] = a_msg[seed_i].read_id;
		if (f_msg->srand[frag_i] == 1)	//+ srand
			f_msg->ex_ref_end[frag_i] = f_msg->ref_end[frag_i] = a_msg[seed_i].offset[aln_i] + seed_len - 1;
		else							//- srand
			f_msg->ex_ref_end[frag_i] = f_msg->ref_end[frag_i] = a_msg[seed_i].offset[aln_i];
		f_msg->ex_read_end[frag_i] = f_msg->read_end[frag_i] = 2*seed_len*(a_msg[seed_i].read_id-1)+seed_len;
		f_msg->per_n[frag_i] = a_msg[seed_i].read_id - a_msg[seed_i].read_id + 1;
	}
	else	//start
	{
		//f_msg->ref_1[frag_i] = a_msg[seed_i].offset[aln_i];
		//f_msg->read_1[frag_i] = a_msg[seed_i].read_id;
		f_msg->frag_num = frag_i + 1;
		f_msg->ex_read_begin[frag_i] = f_msg->read_begin[frag_i] = 2*seed_len*(a_msg[seed_i].read_id-1)+1;
		if (f_msg->srand[frag_i] == 1)	//+ srand
		{
			f_msg->ex_ref_begin[frag_i] = f_msg->ref_begin[frag_i] = a_msg[seed_i].offset[aln_i];
		}
		else							//- srand
		{
			f_msg->ex_ref_begin[frag_i] = f_msg->ref_begin[frag_i] = a_msg[seed_i].offset[aln_i] + seed_len -1;
		}
	}
	return 0;
}

int frag_refresh(frag_msg* f_msg, int f_i, int ref_offset, int read_offset, int FLAG)
{
	if (FLAG == RIGHT)	//right offset
	{
		f_msg->ex_read_end[f_i] += read_offset;
		if (f_msg->srand[f_i] == 1)	//+ srand
			f_msg->ex_ref_end[f_i] += ref_offset;
		else						//- srand
			f_msg->ex_ref_end[f_i] -= ref_offset;
//		f_msg->per_n[f_i] += read_offset / PER_LEN;
	}
	else				//left offset
	{
		f_msg->ex_read_begin[f_i] -= read_offset;
		if (f_msg->srand[f_i] == 1)	//+ srand
			f_msg->ex_ref_begin[f_i] -= ref_offset;
		else
			f_msg->ex_ref_begin[f_i] += ref_offset;
//		f_msg->per_n[f_i] += read_offset / PER_LEN;
	}
	printf("Refresh: %d %d %d %d\n", f_i, ref_offset, read_offset, FLAG);
	return 0;
}

int frag_covered(frag_msg *f_msg, int f_i, int b_f)
{
	f_msg->flag[b_f] = COVERED;
	f_msg->b_f[f_i] = b_f;

	printf("Covered: %d C %d\n", f_i, b_f);
	return 0;
}

/*
 * check long OR short and set ref_len/read_len
 * frag : end   -> return 2
 * long : short -> return 1
 * short : short |
 * long : long  -> return 0
 * short : long -> return -1
 * head : frag  -> return -2
 */
/*
int ls_frag(frag_msg* f_msg, int f_i, int* ref_len, int* read_len)
{
	if (f_i == f_msg->frag_num)	//check if there is un-mapped region before the first frag
	{
		if (f_msg->read_begin[f_msg->c_i[f_i-1]] != 1)	//Bingo!
		{
			*ref_len = f_msg->read_begin[f_msg->c_i[f_i-1]] - 1;
			*read_len = *ref_len;
			return -2;
		}
		return -3;
	}
	else if (f_i == 0)	//last frag
	{
		(*read_len) = f_msg->last_len;
		(*ref_len) = *read_len;
		return 2;
	}

	(*read_len) = f_msg->read_begin[f_msg->c_i[f_i-1]] - f_msg->read_end[f_msg->c_i[f_i]] - 1;
	(*ref_len) = (*read_len);
	
	int ret = 0;

	ret+=(f_msg->per_n[f_msg->c_i[f_i]] < 5 ?-1:0);
	ret+=(f_msg->per_n[f_msg->c_i[f_i-1]] < 5 ?1:0);

	return ret;
}*/

int unsoap2dp_extend(frag_msg *f_msg, int f_i, int FLAG, bntseq_t *bns, int8_t *pac, char *read_seq, int8_t *seq1, int8_t *seq2)
{
	int ret;
	int ref_len, read_len;
	int ref_l_os, read_l_os, ref_r_os, read_r_os;
	int b_f = f_msg->b_f[f_i];
	int ff, m;
	char *seq_s2 = malloc(16385*sizeof(char));

	// XXX Is it necessary to take a 'COVERED frag' as an individaul frag in the next reverse extend?
	//	   Here is YES.
	
	if (FLAG == RIGHT)
	{
		if (b_f == 0)
		{
			read_len = f_msg->last_len;
			ref_len = read_len;
		}
		else
		{
			read_len = f_msg->read_begin[b_f-1] - f_msg->read_end[b_f] - 1;
			ref_len = read_len;
		}
		//XXX pac2fa: left -> OR right ->
		if (f_msg->srand[b_f] == -1)	//- srand
			pac2fa_core(bns, pac, f_msg->chr[b_f], f_msg->ref_end[b_f]-1-ref_len, &ref_len, f_msg->srand[b_f], &ff, seq1);
		else
			pac2fa_core(bns, pac, f_msg->chr[b_f], f_msg->ref_end[b_f], &ref_len, f_msg->srand[b_f], &ff, seq1);
		strncpy(seq_s2, read_seq+f_msg->read_end[b_f], read_len);
		for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
		ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
		frag_refresh(f_msg, f_i, ref_r_os, read_r_os, RIGHT);
	}
	else	//LEFT
	{
		if (b_f == f_msg->frag_num - 1)
		{
			if (f_msg->read_begin[b_f] != 1)	//Un-soap2dp HEAD
			{
				read_len = f_msg->read_begin[b_f] - 1;
				ref_len =read_len;
			}
			else
				return BAD_ALIGN;
		}
		else
		{
			read_len = f_msg->read_begin[b_f] - f_msg->read_end[b_f+1] - 1;
			ref_len = read_len;
		}
		if (f_msg->srand[b_f] == -1)
			pac2fa_core(bns, pac, f_msg->chr[b_f], f_msg->ref_begin[b_f], &ref_len, f_msg->srand[b_f], &ff, seq1);
		else
			pac2fa_core(bns, pac, f_msg->chr[b_f], f_msg->ref_begin[b_f]-ref_len-1, &ref_len, f_msg->srand[b_f], &ff, seq1);
		strncpy(seq_s2, read_seq+(f_msg->read_begin[b_f] - read_len -1), read_len);
		for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
		ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
		frag_refresh(f_msg, f_i, ref_l_os, read_l_os, LEFT);
	}
	return ret;
}

//Extend the frag around with f_i
//Extend PER_LEN in a while cycle, stop when get a BAD_ALIGN
//XXX If it's SHORT, extend the whole frag
//XXX Else, extend PER_LEN in a while cycle, stop when get a BAD_ALIGN
int frag_extend(frag_msg *f_msg, int f_i, int FLAG, bntseq_t *bns, int8_t *pac, char *read_seq, int8_t *seq1, int8_t *seq2, int seed_len)
{
	int ret, i;
	int ref_len, read_len;
	int ref_l_os, read_l_os, ref_r_os, read_r_os;
	int b_f = f_msg->b_f[f_i];
	int ff, m;
	char *seq_s2 = malloc(16385*sizeof(char));

	if (FLAG == RIGHT)
	{
		if (b_f == 0)
		{
			printf("BAD\n");
			return BAD_ALIGN;
		}
		ref_len = read_len = seed_len;
		for (i = 0; i < f_msg->per_n[b_f-1]; i++)
		{
			if (f_msg->srand[f_i] == -1)
				pac2fa_core(bns, pac, f_msg->chr[f_i], f_msg->ex_ref_end[f_i]-1-ref_len, &ref_len, f_msg->srand[f_i], &ff, seq1);
			else
				pac2fa_core(bns, pac, f_msg->chr[f_i], f_msg->ex_ref_end[f_i], &ref_len, f_msg->srand[f_i], &ff, seq1);
			strncpy(seq_s2, read_seq + (f_msg->read_begin[b_f-1] - 1 + i * seed_len), read_len);
			for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
			ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
			frag_refresh(f_msg, f_i, ref_r_os, read_r_os, RIGHT);
			if (ret == BAD_ALIGN)
				return BAD_ALIGN;
		}
		frag_covered(f_msg, f_i, b_f-1);
		/*read_len = f_msg->read_end[b_f-1] - f_msg->read_begin[b_f-1] + 1;
		ref_len = read_len;
		strncpy(seq1, ref_seq->seq.s + f_msg->ex_ref_end[f_i], ref_len);
		strncpy(seq2, read_seq->seq.s + f_msg->ex_read_end[f_i], read_len);
		ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
		frag_refresh(f_msg, f_i, ref_r_os, read_r_os, RIGHT);
		if (ret == GOOD_ALIGN)
			frag_covered(f_msg, f_i, b_f-1);*/
	}
	else	//LEFT
	{
		if (b_f == f_msg->frag_num-1)
			return BAD_ALIGN;
		ref_len = read_len = seed_len;
		for (i = 1; i <= f_msg->per_n[b_f+1]; i++)
		{
			if (f_msg->srand[f_i] == -1)
				pac2fa_core(bns, pac, f_msg->chr[f_i], f_msg->ex_ref_begin[f_i], &ref_len, f_msg->srand[f_i], &ff, seq1);
			else
				pac2fa_core(bns, pac, f_msg->chr[f_i], f_msg->ex_ref_begin[f_i]-1-ref_len, &ref_len, f_msg->srand[f_i], &ff, seq1);
			strncpy(seq_s2, read_seq + (f_msg->read_end[b_f+1] - i * seed_len), read_len);
			for (m = 0; m < read_len; ++m) seq2[m] = nst_nt4_table[(int)seq_s2[m]];
			ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
			frag_refresh(f_msg, f_i, ref_l_os, read_l_os, LEFT);
			if (ret == BAD_ALIGN)
				return BAD_ALIGN;
		}
		frag_covered(f_msg, f_i, b_f+1);
		/*read_len = f_msg->read_end[b_f+1] - f_msg->read_begin[b_f+1] + 1;
		ref_len = read_len;
		strncpy(seq1, ref_seq->seq.s + f_msg->ex_ref_begin[b_f+1]-1, ref_len);
		strncpy(seq2, read_seq->seq.s + f_msg->ex_read_begin[b_f+1]-1, read_len);
		ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
		frag_refresh(f_msg, f_i, ref_r_os, read_r_os, LEFT);
		if (ret == GOOD_ALIGN)
			frag_covered(f_msg, f_i, b_f+1);*/
	}

	return GOOD_ALIGN;
}

int frag_check(bntseq_t *bns, int8_t *pac, char *read_prefix, char *read_seq, frag_msg *f_msg, int seed_len)
{
	int i;
//	fprintf(stdout, "read: %s\n", read_seq);
	fprintf(stdout, "frag:\n");
	for (i = f_msg->frag_num-1; i>= 0; i--)
	{
		fprintf(stdout, "ref %d %d %d read %d %d\n", f_msg->chr[i], f_msg->ref_begin[i], f_msg->ref_end[i], f_msg->read_begin[i], f_msg->read_end[i]);
	}
//	int max_len = 16384 ;
//	char *seq1 = (char*)malloc((max_len+1) * sizeof(char));
//	char *seq2 = (char*)malloc((max_len+1) * sizeof(char));
//
	int max_len = 16384 ;
	int8_t *seq1 = (int8_t*)malloc((max_len+1)*sizeof(int8_t));
	int8_t *seq2 = (int8_t*)malloc((max_len+1)*sizeof(int8_t));
	/*
	 * extend twice
	 */
	//extend right side of every frag
	for (i = f_msg->frag_num-1; i >= 0; i--)
	{
		//if (f_msg->flag[i] == COVERED) continue;
		f_msg->b_f[i] = i;
		while (1)
		{
			if (unsoap2dp_extend(f_msg, i, RIGHT, bns, pac, read_seq, seq1, seq2) == BAD_ALIGN)
				break;
			if (frag_extend(f_msg, i, RIGHT, bns, pac, read_seq, seq1, seq2, seed_len) == BAD_ALIGN)
				break;
		}
	}
	//extend left side of every frag
	for (i = 0; i < f_msg->frag_num; i++)
	{
		//if (f_msg->frag[i] == COVERED) continue;
		f_msg->b_f[i] = i;
		while (1)
		{
			if (unsoap2dp_extend(f_msg, i, LEFT, bns, pac, read_seq, seq1, seq2) == BAD_ALIGN)
				break;
			if (frag_extend(f_msg, i, LEFT, bns, pac, read_seq, seq1, seq2, seed_len) == BAD_ALIGN)
				break;
		}
	}

	for (i = f_msg->frag_num-1; i >= 0; i--)
	{
		if (f_msg->flag[i] == UNCOVERED)
			printf("ref: %d %d  read: %d %d\n", f_msg->ex_ref_begin[i], f_msg->ex_ref_end[i], f_msg->ex_read_begin[i], f_msg->ex_read_end[i]);
	}
	free(seq1);
	free(seq2);
	return 0;
}
