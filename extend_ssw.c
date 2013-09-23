/*
 *	extend_ssw.c
 *	extend_ssw: 
 *	Do Extend-SW on fragments and reference.
 *	Parameter:
 *	frag_start frag_end
 *	ref_start start_end
 *
 *	Created by Yan Gao on 08/26/13.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "ssw.h"
#include "frag_check.h"

#define MAX_MIS 0.05

void ssw_write (s_align* a,	char* ref_seq, char* read_seq, int8_t* table) 
{ 
	fprintf(stdout, "optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t", a->score1, a->score2);
	if (a->ref_begin1 + 1) fprintf(stdout, "target_begin: %d\t", a->ref_begin1 + 1);
	fprintf(stdout, "target_end: %d\t", a->ref_end1 + 1);
	if (a->read_begin1 + 1) fprintf(stdout, "query_begin: %d\t", a->read_begin1 + 1);
	fprintf(stdout, "query_end: %d\n\n", a->read_end1 + 1);
	if (a->cigar) 
	{
		int32_t i, c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
		while (e < a->cigarLen || left > 0) {
			int32_t count = 0;
			int32_t q = qb;
			int32_t p = pb;
			fprintf(stdout, "Target: %8d    ", q + 1);
			for (c = e; c < a->cigarLen; ++c) {
				int32_t letter = 0xf&*(a->cigar + c);
				int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
				int32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if (letter == 1) fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(ref_seq + q));
						++ q;
					}
					++ count;
					if (count == 60) goto step2;
				}
			}
step2:
			fprintf(stdout, "    %d\n                    ", q);
			q = qb;
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				int32_t letter = 0xf&*(a->cigar + c);
				int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
				int32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i){ 
					if (letter == 0) {
						if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)])fprintf(stdout, "|");
						else fprintf(stdout, "*");
						++q;
						++p;
					} else {
						fprintf(stdout, "*");
						if (letter == 1) ++p;
						else ++q;
					}
					++ count;
					if (count == 60) {
						qb = q;
						goto step3;
					}
				}
			}
step3:
			p = pb;
			fprintf(stdout, "\nQuery:  %8d    ", p + 1);
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				int32_t letter = 0xf&*(a->cigar + c);
				int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
				int32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) { 
					if (letter == 2) fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(read_seq + p));
						++p;
					}
					++ count;
					if (count == 60) {
						pb = p;
						left = l - i - 1;
						e = (left == 0) ? (c + 1) : c;
						goto end;
					}
				}
			}
			e = c;
			left = 0;
end:
			fprintf(stdout, "    %d\n\n", p);
		}
	}
}

int result_check(s_align* result, int ref_len, int read_len, int* ref_l_os, int* read_l_os, int* ref_r_os, int* read_r_os)
{
	int ret = GOOD_ALIGN;
	
	if ((double)(result->read_begin1) >= (double)(read_len * MAX_MIS))	//LEFT_BAD
	{
		*ref_r_os = *read_r_os = 0;
		if ((double)(result->read_end1) < (double)(read_len * (1-MAX_MIS)))	//AND RIGHT_BAD
		{
			*ref_l_os = *read_l_os = 0;
		}
		else
		{
			*ref_l_os = ref_len - result->ref_begin1;
			*read_l_os = read_len - result->read_begin1;
		}
		ret = BAD_ALIGN;
	}
	else
	{
		if ((double)(result->read_end1) < (double)(read_len * (1-MAX_MIS)))	//RIGHT_BAD
		{
			*ref_l_os = *read_l_os = 0;
			*ref_r_os = result->ref_end1;
			*read_r_os = result->read_end1;
			ret = BAD_ALIGN;
		}
		else	//GOOD_ALIGN
		{
			*ref_l_os = ref_len - result->ref_begin1 + result->read_begin1;
			*ref_r_os = ref_len + result->ref_end1 - result->read_end1;
			*read_l_os = *read_r_os = read_len;
			ret = GOOD_ALIGN;
		}
	}
	
	printf("ref: %d %d\tread: %d %d\n",result->ref_begin1, result->ref_end1, result->read_begin1, result->read_end1);
	printf("ret: %d\n", ret);
	return ret;
}

//int extend_ssw(char* ref_seq, char* read_seq, int ref_len, int read_len, int* ref_l_os, int* read_l_os, int* ref_r_os, int* read_r_os)
int extend_ssw(int8_t *ref_num, int8_t *read_num, int ref_len, int read_len, int* ref_l_os, int* read_l_os, int* ref_r_os, int* read_r_os)
{
	int i;
	printf("ref: %d\n", ref_len);
	for (i = 0; i < ref_len; i++)
		printf("%d ", ref_num[i]);
	printf("\nread: %d\n", read_len);
	for (i = 0; i < read_len; i++)
		printf("%d ", read_num[i]);
	printf("\n");
	int32_t l, m, k, match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;	// default parameters for genome sequence alignment
	s_profile* profile;
	//int8_t* read_num = (int8_t*)malloc(read_len);
	//int8_t* ref_num  = (int8_t*)malloc(ref_len);	//ref and read seq represented in numbers
	s_align* result;

	/* This table is used to transform nucleotide letters into numbers. */
//	int8_t nt_table[128] = {
//		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
//		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
//		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
//		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
//		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
//	};
	
	// initialize scoring matrix for genome sequences
	//  A  C  G  T	N (or other ambiguous code) 
	//  2 -2 -2 -2 	0	A
	// -2  2 -2 -2 	0	C
	// -2 -2  2 -2 	0	G
	// -2 -2 -2  2 	0	T
	//	0  0  0  0  0	N (or other ambiguous code)	
	int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
	for (l = k = 0; l < 4; ++l) {
		for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base: no penalty
	}
	for (m = 0; m < 5; ++m) mat[k++] = 0;

	//for (m = 0; m < read_len; ++m) read_num[m] = nt_table[(int)read_seq[m]];
	profile = ssw_init(read_num, read_len, mat, 5, 2);
	//for (m = 0; m < ref_len; ++m) ref_num[m] = nt_table[(int)ref_seq[m]];

	result = ssw_align (profile, ref_num, ref_len, gap_open, gap_extension, 1, 0, 0, 15);	
	//ssw_write(result, ref_seq, read_seq, nt_table);

	free(mat);
	//free(ref_num);
	//free(read_num);
	
	//check if result represent a good align or not
	/*set offset of every side
	 *	check if result represent a good align or not
	 *	Good : return 0	(align len >= 90%)
	 *	Bad  : return 1
	 */
	return(result_check(result, ref_len, read_len, ref_l_os, read_l_os, ref_r_os, read_r_os));
}
