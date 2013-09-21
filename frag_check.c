/*
 * frag_check.c
 *
 * function: transform raw frag.msg to detailed frag.msg
 * method:   extend-ssw
 *
 * Created by Yan Gao on 08/29/2013
 */

#include "frag_check.h"

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include "ssw.h"
#include "kseq.h"

#define PER_LEN 100

KSEQ_INIT(gzFile, gzread)

/*
 * parameter: frag.msg
 */

int extend_ssw(char* ref_seq, char* read_seq, int ref_len, int read_len, int* ref_l_os, int* read_l_os, int* ref_r_os, int* read_r_os);

frag_msg* frag_init_msg(int frag_max)
{
	frag_msg *f_msg = (frag_msg*)malloc(sizeof(frag_msg));

	f_msg->frag_max = frag_max; 
	f_msg->frag_num = 0;
	f_msg->last_len = PER_LEN;	//XXX
	f_msg->read_num = 500;		//XXX
	f_msg->chr= (int *)malloc(frag_max * sizeof(int));
	f_msg->srand = (int *)malloc(frag_max * sizeof(int));
	f_msg->read_1 = (int *)malloc(frag_max * sizeof(int));
	f_msg->read_2 = (int *)malloc(frag_max * sizeof(int));
	f_msg->ref_1 = (int *)malloc(frag_max * sizeof(int));
	f_msg->ref_2 = (int *)malloc(frag_max * sizeof(int));

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
	free(f_msg->read_1);
	free(f_msg->read_2);
	free(f_msg->ref_1);
	free(f_msg->ref_2);
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

int frag_refresh(frag_msg* f_msg, int f_i, int ref_offset, int read_offset, int FLAG)
{
	if (FLAG == RIGHT)	//right offset
	{
		f_msg->ex_read_end[f_i] += read_offset;
		f_msg->ex_ref_end[f_i] += ref_offset;
//		f_msg->per_n[f_i] += read_offset / PER_LEN;
	}
	else				//left offset
	{
		f_msg->ex_read_begin[f_i] -= read_offset;
		f_msg->ex_ref_begin[f_i] -= ref_offset;
//		f_msg->per_n[f_i] += read_offset / PER_LEN;
	}
	printf("Refresh: %d %d %d %d\n", f_i, ref_offset, read_offset, FLAG);
	return 0;
}
/*
int frag_refresh(frag_msg* f_msg, int f_i, int ref_offset, int read_offset, int FLAG)
{
	if (FLAG == RIGHT)	//right offset
	{
		f_msg->read_end[f_msg->c_i[f_i]] += read_offset;
		f_msg->ref_end[f_msg->c_i[f_i]] += ref_offset;
		f_msg->per_n[f_msg->c_i[f_i]] += read_offset / PER_LEN;
	}
	else				//left offset
	{
		f_msg->read_begin[f_msg->c_i[f_i]] -= read_offset;
		f_msg->ref_begin[f_msg->c_i[f_i]] -= ref_offset;
		f_msg->per_n[f_msg->c_i[f_i]] += read_offset / PER_LEN;
	}
	printf("Refresh: %d %d %d %d\n", f_msg->c_i[f_i], ref_offset, read_offset, FLAG);
	return 0;
}
*/
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

int unsoap2dp_extend(frag_msg* f_msg, int f_i, int FLAG, kseq_t* ref_seq, kseq_t* read_seq, char* seq1, char* seq2)
{
	int ret;
	int ref_len, read_len;
	int ref_l_os, read_l_os, ref_r_os, read_r_os;
	int b_f = f_msg->b_f[f_i];

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
		strncpy(seq1, ref_seq->seq.s+f_msg->ref_end[b_f], ref_len);
		strncpy(seq2, read_seq->seq.s+f_msg->read_end[b_f], read_len);
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
		strncpy(seq1, ref_seq->seq.s+(f_msg->ref_begin[b_f] - ref_len - 1), ref_len);
		strncpy(seq2, read_seq->seq.s+(f_msg->read_begin[b_f] - read_len -1), read_len);
		ret = extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os);
		frag_refresh(f_msg, f_i, ref_l_os, read_l_os, LEFT);
	}
	return ret;
}

//Extend the frag around with f_i
//Extend PER_LEN in a while cycle, stop when get a BAD_ALIGN
//XXX If it's SHORT, extend the whole frag
//XXX Else, extend PER_LEN in a while cycle, stop when get a BAD_ALIGN
int frag_extend(frag_msg* f_msg, int f_i, int FLAG, kseq_t* ref_seq, kseq_t* read_seq, char* seq1, char* seq2)
{
	int ret, i;
	int ref_len, read_len;
	int ref_l_os, read_l_os, ref_r_os, read_r_os;
	int b_f = f_msg->b_f[f_i];

	if (FLAG == RIGHT)
	{
		if (b_f == 0)
		{
			printf("BAD\n");
			return BAD_ALIGN;
		}
		ref_len = read_len = PER_LEN;
		for (i = 0; i < f_msg->per_n[b_f-1]; i++)
		{
			strncpy(seq1, ref_seq->seq.s + f_msg->ex_ref_end[f_i], ref_len);
			strncpy(seq2, read_seq->seq.s + (f_msg->read_begin[b_f-1] - 1 + i * PER_LEN), read_len);
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
		ref_len = read_len = PER_LEN;
		for (i = 1; i <= f_msg->per_n[b_f+1]; i++)
		{
			strncpy(seq1, ref_seq->seq.s + f_msg->ex_ref_begin[f_i] - 1 - ref_len, ref_len);
			strncpy(seq2, read_seq->seq.s + (f_msg->read_end[b_f+1] - i * PER_LEN), read_len);
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

int frag_check(char *frag_f)
{
	FILE *fp = fopen(frag_f, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "ERROR, can't open %s.", frag_f);
		return 1;
	}
	frag_msg *f_msg = frag_init_msg(MAX_FRAG);
	char line[MAX_LINE_LEN+1];
	int n = 0;	//frag num
	//
	//get ref and contig seq
	char ref[100] = "../data/chr1_hg18.fa";
	char read[100]= "../data/venter_chr1_100000_0.01_het.fq";

	//get frag msg
	while (fgets(line, MAX_LINE_LEN, fp) != NULL)
	{
		sscanf(line, "%d\t%d\t%d\t%d\t%d\t%d", f_msg->chr+n, f_msg->srand+n, f_msg->ref_2+n, f_msg->read_2+n, f_msg->ref_1+n, f_msg->read_1+n);
		f_msg->frag_num++;
		f_msg->ex_read_begin[n] = f_msg->read_begin[n] = 2 * PER_LEN * (f_msg->read_1[n] - 1) + 1;
		f_msg->ex_read_end[n] = f_msg->read_end[n] = 2 * PER_LEN * (f_msg->read_2[n] - 1) + PER_LEN;
		f_msg->ex_ref_begin[n] = f_msg->ref_begin[n] =  f_msg->ref_1[n];
		f_msg->ex_ref_end[n] = f_msg->ref_end[n] = f_msg->ref_2[n] + PER_LEN - 1;
		f_msg->per_n[n] = f_msg->read_2[n] - f_msg->read_1[n] + 1;
		printf("%d %d %d %d %d %d %d\n", f_msg->frag_num, f_msg->chr[n], f_msg->srand[n], f_msg->ref_2[n], f_msg->read_2[n], f_msg->ref_1[n], f_msg->read_1[n]);
		n++;
		//		if (f_msg->frag_num == f_msg->frag_max)
		//		frag_realloc(f_msg);	//double mem
	}
	fclose(fp);


	int max_len = 16384 ;
	char *seq1 = (char*)malloc((max_len+1) * sizeof(char));
	char *seq2 = (char*)malloc((max_len+1) * sizeof(char));
	gzFile reffp, readfp;
	kseq_t *ref_seq, *read_seq;

	int chr=0, ctg=0;

	reffp = gzopen(ref, "r");
	readfp = gzopen(read, "r");
	ref_seq = kseq_init(reffp);
	read_seq = kseq_init(readfp);

	while (kseq_read(ref_seq) > 0)
	{
		chr++;	//get one chromosome
		break;
	}
	

	while(kseq_read(read_seq) > 0)
	{
		ctg++;
		if (ctg == 1620)
			break;
	}

	//strncpy(seq2, read_seq->seq.s+20600, read_len);
	
	/*
	int i;
	int ret, ref_l_os=0, read_l_os=0, ref_r_os=0, read_r_os=0;	//left_os right_os : offset of align result 

	for (i = f_msg->frag_num; i >= 0; i--)
	{
		//long & short || long & long, short & short
		
		printf("%d ", i);
		ret = ls_frag(f_msg, i, &ref_len, &read_len);
		printf("%d %d %d\n", ret, ref_len, read_len);

		if (ret == -3)	//
			continue;
		if (ret >= 0)	// long : short	| both
		{
			//extend un-soap2dp PER_LEN
			strncpy(seq1, ref_seq->seq.s + f_msg->ref_end[f_msg->c_i[i]], ref_len);
			strncpy(seq2, read_seq->seq.s + f_msg->read_end[f_msg->c_i[i]], read_len);

			if (extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os) == GOOD_ALIGN)	//extend the frag behind
			{
				frag_refresh(f_msg, i, ref_r_os, read_r_os, RIGHT);	//refresh frag with un-soap2dp region
				if (ret == 1)			//next frag is SHORT, extend all the frag
				{
					read_len = f_msg->read_end[f_msg->c_i[i-1]] - f_msg->read_begin[f_msg->c_i[i-1]] + 1;
					ref_len = read_len;
					strncpy(seq1, ref_seq->seq.s + f_msg->ref_end[f_msg->c_i[i]], ref_len);
					strncpy(seq2, read_seq->seq.s + f_msg->read_end[f_msg->c_i[i]], read_len);
					if (extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os) == GOOD_ALIGN)	
					//refresh frag with all the next frag, mark this short frag as 'covered'
					{
						frag_refresh(f_msg, i, ref_r_os, read_r_os, RIGHT);
						frag_covered(f_msg, i-1);
					}
					else	//
					{
						frag_refresh(f_msg, i, ref_r_os, read_r_os, RIGHT);	
					}
				}
				else if (ret == 0)		//next frag is LONG, extend one PER_LEN
				{
					ref_len = read_len = PER_LEN;
					strncpy(seq1, ref_seq->seq.s + f_msg->ref_end[f_msg->c_i[i]], ref_len);
					strncpy(seq2, read_seq->seq.s + f_msg->read_end[f_msg->c_i[i]], read_len);

					while (extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os) == GOOD_ALIGN)
					{
						frag_refresh(f_msg, i, ref_r_os, read_r_os, RIGHT);
						strncpy(seq1, ref_seq->seq.s + f_msg->ref_end[f_msg->c_i[i]], ref_len);
						strncpy(seq2, read_seq->seq.s + f_msg->read_end[f_msg->c_i[i]], read_len);
					}
					frag_refresh(f_msg, i, ref_r_os, read_r_os, RIGHT);
				}
				//else					//ret == 2, no next frag
				//{
				//	frag_refresh(f_msg, i, ref_r_os, read_r_os, RIGHT);
				//}
			}
			else																				//Bad align
				frag_refresh(f_msg, i, ref_r_os, read_r_os, RIGHT);
		}
		if (ret <= 0)	// short : long | both
		{
			strncpy(seq1, ref_seq->seq.s + f_msg->ref_begin[f_msg->c_i[i-1]] - ref_len - 1, ref_len);
			strncpy(seq2, read_seq->seq.s + f_msg->read_begin[f_msg->c_i[i-1]] - read_len - 1, read_len);

			if (extend_ssw(seq1, seq2, ref_len , read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os) == GOOD_ALIGN)	//extend the frag before
			{
				frag_refresh(f_msg, i-1, ref_l_os, read_l_os, LEFT);	//refresh frag with un-soap2dp region
				if (ret == -1)	//previous frag is SHORT, extend all the frag
				{
					read_len = f_msg->read_end[f_msg->c_i[i]] - f_msg->read_begin[f_msg->c_i[i]] + 1;
					ref_len = read_len;
					strncpy(seq1, ref_seq->seq.s + f_msg->ref_begin[f_msg->c_i[i-1]] - ref_len -1, ref_len);
					strncpy(seq2, read_seq->seq.s + f_msg->read_begin[f_msg->c_i[i]], read_len);
					if (extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os) == GOOD_ALIGN)
					//refresh frag with the whole previous frag, mark this short frag as 'covered'
					{
						frag_refresh(f_msg, i-1, ref_l_os, read_l_os, LEFT);
						frag_covered(f_msg, i);

					}
					else
						frag_refresh(f_msg, i-1, ref_l_os, read_l_os, LEFT);
				}
				else if (ret == 0)			//previous frag is LONG, extend one PER_LEN
				{
					ref_len = read_len = PER_LEN;
					strncpy(seq1, ref_seq->seq.s + f_msg->ref_begin[f_msg->c_i[i-1]] - ref_len - 1, ref_len);
					strncpy(seq2, read_seq->seq.s + f_msg->read_begin[f_msg->c_i[i-1]] - read_len -1, read_len);

					while (extend_ssw(seq1, seq2, ref_len, read_len, &ref_l_os, &read_l_os, &ref_r_os, &read_r_os) == GOOD_ALIGN)
					{
						frag_refresh(f_msg, i-1, ref_l_os, read_l_os, LEFT);
						strncpy(seq1, ref_seq->seq.s + f_msg->ref_begin[f_msg->c_i[i-1]] - ref_len - 1, ref_len);
						strncpy(seq2, read_seq->seq.s + f_msg->read_begin[f_msg->c_i[i-1]] - read_len -1, read_len);
					}
					frag_refresh(f_msg, i-1, ref_l_os, read_l_os, LEFT);
				}
				//else						//ret == -2, no previous frag
				//{
				//	frag_refresh(f_msg, i, ref_l_os, read_l_os, LEFT);
				//}
			}
			else															//Bad align
				frag_refresh(f_msg, i-1, ref_l_os, read_l_os, LEFT);
		}
	}
	*/

	/*
	 * extend twice
	 */
	//extend right side of every frag
	//
	int i;
	for (i = f_msg->frag_num-1; i >= 0; i--)
	{
		//if (f_msg->flag[i] == COVERED)
		//	continue;
		f_msg->b_f[i] = i;
		while (1)
		{
			if (unsoap2dp_extend(f_msg, i, RIGHT, ref_seq, read_seq, seq1, seq2) == BAD_ALIGN)
				break;
			if (frag_extend(f_msg, i, RIGHT, ref_seq, read_seq, seq1, seq2) == BAD_ALIGN)
				break;
		}
	}
	
	//extend left side of every frag
	for (i = 0; i < f_msg->frag_num; i++)
	{
		//if (f_msg->frag[i] == COVERED)
		//	continue;
		f_msg->b_f[i] = i;
		while (1)
		{
			if (unsoap2dp_extend(f_msg, i, LEFT, ref_seq, read_seq, seq1, seq2) == BAD_ALIGN)
				break;
			if (frag_extend(f_msg, i, LEFT, ref_seq, read_seq, seq1, seq2) == BAD_ALIGN)
				break;
		}
	}

	for (i = f_msg->frag_num-1; i >= 0; i--)
	{
		if (f_msg->flag[i] == UNCOVERED)
			printf("ref: %d %d  read: %d %d\n", f_msg->ex_ref_begin[i], f_msg->ex_ref_end[i], f_msg->ex_read_begin[i], f_msg->ex_read_end[i]);
	}

//	frag_free(f_msg);
	free(seq1);
	free(seq2);
	gzclose(reffp);
	gzclose(readfp);

	return 0;
}
