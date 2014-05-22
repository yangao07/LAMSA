#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include "bam.h"
#include "sam_header.h"
#include "sam.h"
#include "kstring.h"
#include "khash.h"
#include "kseq.h"
#include "../lsat_aln.h"


KSTREAM_INIT(gzFile, gzread, 16384)
KHASH_SET_INIT_STR(rg)

samfile_t *samopen(const char *fn, const char *mode)
{
	samfile_t *fp;
	fp = (samfile_t*)calloc(1, sizeof(samfile_t));
	if (strchr(mode, 'r')) { // read

		{ // text
			fp->x.tamr = sam_open(fn);
			if (fp->x.tamr == 0) goto open_err_ret;
			fp->header = sam_header_read(fp->x.tamr);
			/*if (fp->header->n_targets == 0) { // no @SQ fields
				if (fp->header->n_targets == 0 && bam_verbose >= 1)
					fprintf(stderr, "[samopen] no @SQ lines in the header.\n");
			} else if (bam_verbose >= 2) fprintf(stderr, "[samopen] SAM header is present: %d sequences.\n", fp->header->n_targets);*/
		}
	} 
	return fp;

open_err_ret:
	free(fp);
	return 0;
}


void samclose(samfile_t *fp)
{
	if (fp == 0) return;
	if (fp->header) bam_header_destroy(fp->header);
	sam_close(fp->x.tamr);
	free(fp);
}

//XXX start
/*int sam_view_main(int argc, char *argv[])
{
	int ret = 0;
	samfile_t *in = 0;

	if ((in = samopen(argv[optind], "r")) == 0) {
		fprintf(stderr, "[main_samview] fail to open \"%s\" for reading.\n", argv[optind]);
		ret = 1;
		goto view_end;
	}
	if (in->header == 0) {
		fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
		ret = 1;
		goto view_end;
	}

	//if (n_threads > 1) samthreads(out, n_threads, 256); 
	//if (is_header_only) goto view_end; // no need to print alignments

	{ // convert/print the entire file
		bam1_t *b = bam_init1();
		int r;
		while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
			;
			}
		if (r < -1) {
			fprintf(stderr, "[main_samview] truncated file.\n");
			ret = 1;
		}
		bam_destroy1(b);
	} 
view_end:
	// close files, free and return
	samclose(in);
	return ret;
}*/

/*
int main(int argc, char *argv[])
{
	int i;
	int ret;
	samfile_t *in = 0;
	if ((in = samopen(argv[optind], "r")) == 0) {
		fprintf(stderr, "[main_samview] fail to open \"%s\" for reading.\n", argv[optind]);
		ret = 1;
		goto view_end;
	}
	if (in->header == 0) {
		fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
		ret = 1;
		goto view_end;
	}
	
	sam_msg *m_msg = (sam_msg*)calloc(1, sizeof(sam_msg));
	m_msg->sam_n = 0; m_msg->sam_m = 1;
	m_msg->sam = (sam_t*)malloc(sizeof(sam_t));
	m_msg->sam->cigar_s = (kstring_t*)calloc(1, sizeof(kstring_t));

	int r;
	int xa_limit = 100;
	
	while ((r = sam_read1(in->x.tamr, in->header, m_msg, xa_limit)) >= 0) { // read one alignment from `in'
		printf("sam:\n");
		for (i = 0; i < r; ++i)
			printf("\t%d\t%d\t%lld\t%s\n", m_msg->sam[i].nsrand, m_msg->sam[i].chr, (long long)m_msg->sam[i].offset, m_msg->sam[i].cigar_s->s);
	}
	if (r < -1) {
		fprintf(stderr, "[main_samview] truncated file.\n");
		ret = 1;
	}

	for (i = 0; i < in->header->n_targets; ++i)
		fprintf(stdout, "t %d:\t%s\n", i+1, in->header->target_name[i]);

view_end:
	// close files, free and return
	for (i = 0; i < m_msg->sam_m; ++i)
	{	
		if(m_msg->sam[i].cigar_s->s) free(m_msg->sam[i].cigar_s->s);
		free(m_msg->sam[i].cigar_s);
	}
	free(m_msg->sam); free(m_msg);
	samclose(in);
	return ret;
}*/
//XXX end
