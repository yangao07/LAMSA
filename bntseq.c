/*
 * bntseq.c
 * Transform fa files to binary files.
 * Output: prefix.lsat.pac. 
 * Use '00'/'01'/'10'/'11' to represent 'A'/'G'/'C'/'T'.
 * For short 'N's(<= 10bp), replaced by 'G'(same as 'soap2-dp').
 * For long 'N's(> 10bp), just throw them away when mapping.
 * 0-base coordinate. 
 */

/*
 * Created by Y Gao on 2013/9/16.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bntseq.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

char n_char[6] = {'A', 'C', 'G', 'T', 'N' };

void bns_dump(const bntseq_t *bns, const char *prefix)
{
	char str[1024];
	FILE *fp;
	int i;
	{ // dump .ann
		strcpy(str, prefix); strcat(str, ".lsat.ann");
		fp = fopen(str, "w");
		fprintf(fp, "%lld %d\n", (long long)bns->l_pac, bns->n_seqs);
		for (i = 0; i != bns->n_seqs; ++i) 
		{
			bntann1_t *p = bns->anns + i;
			fprintf(fp, "%d %s", p->gi, p->name);
			if (p->anno[0]) fprintf(fp, " %s\n", p->anno);
			else fprintf(fp, "\n");
			fprintf(fp, "%lld %d %d %d\n", (long long)p->offset, p->len, p->ambs_offset, p->n_ambs);
		}
		fclose(fp);
	}
	{ // dump .amb
		strcpy(str, prefix); strcat(str, ".lsat.amb");
		fp = fopen(str, "w");
		fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->n_holes);
		for (i = 0; i != bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fprintf(fp, "%lld %d %c\n", (long long)p->offset, p->len, p->amb);
		}
		fclose(fp);
	}
}

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else
	{
		int i;
		if (bns->fp_pac)	fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i != bns->n_seqs; ++i)
		{
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

void bns_fa2bnt(gzFile fp_fa, const char *prefix)
{
	kseq_t *seq;
	char name[1024];
	bntseq_t *bns;
	bntamb1_t *q;
	int l_buf;
	char buf[0x10000];
	int32_t m_seqs, m_holes, l, i;
	FILE *fp;

	// initialization
	seq = kseq_init(fp_fa);
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	m_seqs = m_holes = 8;
	bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
	bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
	q = bns->ambs;
	l_buf = 0;
	strcpy(name, prefix); strcat(name, ".lsat.pac");
	fp = fopen(name, "wb");
	memset(buf, 0, 0x10000);
	// read sequences
	while ((l = kseq_read(seq)) >= 0) 
	{
		bntann1_t *p;
		int lasts;
		if (bns->n_seqs == m_seqs) 
		{
			m_seqs <<= 1;
			bns->anns = (bntann1_t*)realloc(bns->anns, m_seqs * sizeof(bntann1_t));
		}
		p = bns->anns + bns->n_seqs;
		p->name = strdup((char*)seq->name.s);
		p->anno = seq->comment.s? strdup((char*)seq->comment.s) : strdup("(null)");
		p->gi = 0; p->len = l;
		p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
		p->n_ambs = 0;
		p->ambs_offset = (bns->n_seqs == 0)? 0 : (p-1)->ambs_offset + (p-1)->n_ambs;
		for (i = 0, lasts = 0; i < l; ++i) 
		{
			int c = nst_nt4_table[(int)seq->seq.s[i]];
			if (c >= 4) 
			{ // N
				if (lasts == seq->seq.s[i])  // contiguous N
					++q->len;
				else 
				{
					if (bns->n_holes == m_holes) 
					{
						m_holes <<= 1;
						bns->ambs = (bntamb1_t*)realloc(bns->ambs, m_holes * sizeof(bntamb1_t));
					}
					q = bns->ambs + bns->n_holes;
					q->len = 1;
					q->offset = p->offset + i;
					q->amb = seq->seq.s[i];
					++p->n_ambs;
					++bns->n_holes;
				}
			}
			lasts = seq->seq.s[i];
			{ // fill buffer
				if (c >= 4) c = 2;	//'G'
				if (l_buf == 0x40000) 
				{
					fwrite(buf, 1, 0x10000, fp);
					memset(buf, 0, 0x10000);
					l_buf = 0;
				}
				buf[l_buf>>2] |= c << ((3 - (l_buf&3)) << 1);
				++l_buf;
			}
		}
		++bns->n_seqs;
		bns->l_pac += seq->seq.l;
	}
	{ // finalize .pac file
		char ct;
		fwrite(buf, 1, (l_buf>>2) + ((l_buf&3) == 0? 0 : 1), fp);
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % 4 == 0) 
		{
			ct = 0;
			fwrite(&ct, 1, 1, fp);
		}
		ct = bns->l_pac % 4;
		fwrite(&ct, 1, 1, fp);
		// close .pac file
		fclose(fp);
	}
	bns_dump(bns, prefix);
	bns_destroy(bns);
	kseq_destroy(seq);
}

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[1024];
	FILE *fp;
	bntseq_t *bns;
	long long xx;
	int i;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = fopen(ann_filename, "r");
		fscanf(fp, "%lld%d", &xx, &bns->n_seqs);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments 
			while ((c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			fscanf(fp, "%lld%d%d%d", &xx, &p->len, &p->ambs_offset, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = fopen(amb_filename, "r");
		fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		l_pac = xx;
		//XXX//xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = fopen(pac_filename, "rb");
	}
	return bns;
}

bntseq_t *bns_restore(const char *prefix)
{  
	char ann_filename[1024], amb_filename[1024], pac_filename[1024];
	strcat(strcpy(ann_filename, prefix), ".lsat.ann");
	strcat(strcpy(amb_filename, prefix), ".lsat.amb");
	strcat(strcpy(pac_filename, prefix), ".lsat.pac");
	return bns_restore_core(ann_filename, amb_filename, pac_filename);
}

int32_t n_recover(const bntseq_t *bns, const int32_t seq_n, const int64_t pac_coor, const int64_t len, int8_t *seq, const int srand)
{
	int i;
	int32_t offset = bns->anns[seq_n-1].ambs_offset;

	// binary search for most left and most right ambs
	int left, mid, right, l_amb, r_amb;
	left = 0; right = bns->anns[seq_n-1].n_ambs;
	l_amb = r_amb = -1;
	while (left < right)	//search for most left ambs
	{
		mid = (left + right) >> 1;
		if (bns->ambs[mid+offset].offset+bns->ambs[mid+offset].len <= pac_coor)
			left = mid + 1;
		else if (bns->ambs[mid+offset].offset >= pac_coor + len)
			right = mid;
		else	//overlap
		{
			if (mid == 0 || bns->ambs[mid+offset-1].offset+bns->ambs[mid+offset-1].len <= pac_coor)	//most left
			{
				r_amb = l_amb = mid;
				break;
			}
			else
				right = mid;
		}
	}
	if (l_amb == -1)	//no 'N'
		return 0;
	
	left = l_amb + 1; right = bns->anns[seq_n-1].n_ambs;
	while (left < right)
	{
		mid = (left + right) >> 1;
		if (bns->ambs[mid+offset].offset >= pac_coor + len)
			right = mid;
		else 
		{
			if (mid == right-1 || bns->ambs[mid+1+offset].offset >= pac_coor + len)	//right most
			{
				r_amb = mid;
				break;
			}
			else
				left = mid + 1;
		}
	}
	for (i = l_amb; i <= r_amb; i++)
	{
		int64_t s, e, j;
		s = pac_coor > bns->ambs[i+offset].offset ? pac_coor : bns->ambs[i+offset].offset;
		e = pac_coor + len > (bns->ambs[i+offset].offset + bns->ambs[i+offset].len) ? 
				bns->ambs[i+offset].offset + bns->ambs[i+offset].len : pac_coor+len;
		printf("s: %lld e: %lld\n", s, e);
		if (srand == -1)	//rev
		{
			for (j = len-e; j < len-s; j++)
				seq[j+pac_coor] <<= 2;
		}
		else
		{
			for (j = s; j < e; j++)
				seq[j-pac_coor] <<= 1;	//recover to N
		}
	}
	return 2;
}

#define __rpac(pac, l, i) (pac[(l-i-1)>>2] >> (~(l-i-1)&3)*2 & 0x3)
//convert binary file to ACGT sequence
int pac2fa_core(const bntseq_t *bns, const int8_t *pac, const int32_t seq_n, const int64_t start/*0-base*/, int32_t *len, const int srand, int *FLAG, int8_t *seq)
{
	int64_t pac_coor;
	int i,k;
	
	*FLAG = 0;

	if (start > bns->anns[seq_n-1].len)
	{
		fprintf(stderr, "[bntseq] ERROR: Coor is longger than sequence lenth.(%lld > %d)\n", (long long)start, bns->anns[seq_n-1].len);
		exit(1);
	}
	pac_coor = bns->anns[seq_n-1].offset + start;
	if (start + *len > bns->anns[seq_n-1].len)	//candidate seq is out of range.
	{
		*len = bns->anns[seq_n-1].len - start;
		*FLAG |= 1;	
	}

	if (srand == -1)	//rev
	{
		for (i = *len-1, k = pac_coor; i >= 0; i--, k++)
			seq[i] = 3-(pac[k>>2] >> ((~k&3) << 1) & 0x3);
	}
	else
	{
		for (i = 0, k = pac_coor; i < *len; i++, k++)
			seq[i] = (pac[k>>2] >> ((~k&3) << 1) & 0x3);
	}

	*FLAG |= n_recover(bns, seq_n, pac_coor, *len, seq, srand);	//0 : no 'N', or 2.
	return 0;
}

//int main(int argc, char *argv[])
int fa2pac(int argc, char *argv[])
{
	gzFile fp;

	if (argc < 2)
	{
		fprintf(stderr, "Usage: fa2pac <in.fa> [<out.prefix>]\n");
		return 1;
	}
	fp = gzopen(argv[1], "r");
	bns_fa2bnt(fp, (argc < 3)? argv[1] : argv[2]);
	gzclose(fp);

	//test for pac2fa
	bntseq_t *bns;
	bns = bns_restore(argv[1]);
	int8_t *pac, *seq;

	pac = (int8_t*)calloc(bns->l_pac/4 + 1, 1);
	fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);

	int32_t len = 260;
	int flag=0,i;
	seq = (int8_t*)calloc(len, 1);
	pac2fa_core(bns, pac, 1, 0, &len, 1, &flag, seq);

	printf("len:  %d\nflag: %d\n", len, flag);
	for (i = 0; i < len; i++)
	{
		printf("%c", n_char[seq[i]]);
		if ((i+1) % 50 == 0)
			printf("\n");
	}
	printf("\n");
	free(pac);
	free(seq);
	bns_destroy(bns);

	return 0;
}
