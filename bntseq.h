#ifndef BNTSEQ_H
#define BNTSEQ_H

#include <stdint.h>
#include <zlib.h>

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t ambs_offset;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
//	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;


extern char nst_nt4_table[256];
//int fa2pac(int argc, char *argv[]);

void bns_fa2bnt(gzFile fp_fa, const char *prefix);
bntseq_t *bns_restore(const char *prefix);
void bns_destroy(bntseq_t *bns);
int pac2fa_core(const bntseq_t *bns, const uint8_t *pac, const int32_t seq_n, const int64_t start/*0-base*/, int32_t *len, const int rev, int *FLAG, uint8_t *seq);
#endif
