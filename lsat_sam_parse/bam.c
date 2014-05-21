#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include "bam.h"
#include "kstring.h"
#include "sam_header.h"
#include "bam_endian.h"

int bam_is_be = 0, bam_verbose = 2, bam_no_B = 0;
char *bam_flag2char_table = "pPuUrR12sfd\0\0\0\0\0";

/**************************
 * CIGAR related routines *
 **************************/

//XXX s

uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar)
{
	int k, end = c->pos;
	for (k = 0; k < c->n_cigar; ++k) {
		int op  = bam_cigar_op(cigar[k]);
		int len = bam_cigar_oplen(cigar[k]);
		if (op == BAM_CBACK) { // move backward
			int l, u, v;
			if (k == c->n_cigar - 1) break; // skip trailing 'B'
			for (l = k - 1, u = v = 0; l >= 0; --l) {
				int op1  = bam_cigar_op(cigar[l]);
				int len1 = bam_cigar_oplen(cigar[l]);
				if (bam_cigar_type(op1)&1) { // consume query
					if (u + len1 >= len) { // stop
						if (bam_cigar_type(op1)&2) v += len - u;
						break;
					} else u += len1;
				}
				if (bam_cigar_type(op1)&2) v += len1;
			}
			end = l < 0? c->pos : end - v;
		} else if (bam_cigar_type(op)&2) end += bam_cigar_oplen(cigar[k]);
	}
	return end;
}

int32_t bam_cigar2qlen(const bam1_core_t *c, const uint32_t *cigar)
{
	uint32_t k;
	int32_t l = 0;
	for (k = 0; k < c->n_cigar; ++k)
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			l += bam_cigar_oplen(cigar[k]);
	return l;
}

/********************
 * BAM I/O routines *
 ********************/

bam_header_t *bam_header_init()
{
	bam_is_be = bam_is_big_endian();
	return (bam_header_t*)calloc(1, sizeof(bam_header_t));
}

void bam_header_destroy(bam_header_t *header)
{
	int32_t i;
	extern void bam_destroy_header_hash(bam_header_t *header);
	if (header == 0) return;
	if (header->target_name) {
		for (i = 0; i < header->n_targets; ++i)
			free(header->target_name[i]);
		free(header->target_name);
		free(header->target_len);
	}
	free(header->text);
	if (header->dict) sam_header_free(header->dict);
	if (header->rg2lib) sam_tbl_destroy(header->rg2lib);
	bam_destroy_header_hash(header);
	free(header);
}

// XXX e
