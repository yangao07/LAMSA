#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#ifdef _WIN32
#include <fcntl.h>
#endif
#include "kstring.h"
#include "bam.h"
#include "sam_header.h"
#include "kseq.h"
#include "khash.h"
#include "../lsat_aln.h"

KSTREAM_INIT(gzFile, gzread, 16384)
KHASH_MAP_INIT_STR(ref, uint64_t)

void bam_init_header_hash(bam_header_t *header);
void bam_destroy_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

// XXX s
unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};
// XXX e

unsigned short bam_char2flag_table[256] = {
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,BAM_FREAD1,BAM_FREAD2,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	BAM_FPROPER_PAIR,0,BAM_FMREVERSE,0, 0,BAM_FMUNMAP,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, BAM_FDUP,0,BAM_FQCFAIL,0, 0,0,0,0, 0,0,0,0,
	BAM_FPAIRED,0,BAM_FREVERSE,BAM_FSECONDARY, 0,BAM_FUNMAP,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
};

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

struct __tamFile_t {
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	uint64_t n_lines;
	int is_first;
};

//XXX s
static bam_header_t *hash2header(const kh_ref_t *hash)
{
	bam_header_t *header;
	khiter_t k;
	header = bam_header_init();
	header->n_targets = kh_size(hash);
	header->target_name = (char**)calloc(kh_size(hash), sizeof(char*));
	header->target_len = (uint32_t*)calloc(kh_size(hash), 4);
	for (k = kh_begin(hash); k != kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			int i = (int)kh_value(hash, k);
			header->target_name[i] = (char*)kh_key(hash, k);
			header->target_len[i] = kh_value(hash, k)>>32;
		}
	}
	bam_init_header_hash(header);
	return header;
}

bam_header_t *sam_header_read2(const char *fn)
{
	bam_header_t *header;
	int c, dret, ret, error = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	kh_ref_t *hash;
	khiter_t k;
	if (fn == 0) return 0;
	fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	if (fp == 0) return 0;
	hash = kh_init(ref);
	ks = ks_init(fp);
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	while (ks_getuntil(ks, 0, str, &dret) > 0) {
		char *s = strdup(str->s);
		int len, i;
		i = kh_size(hash);
		ks_getuntil(ks, 0, str, &dret);
		len = atoi(str->s);
		k = kh_put(ref, hash, s, &ret);
		if (ret == 0) {
			fprintf(stderr, "[sam_header_read2] duplicated sequence name: %s\n", s);
			error = 1;
		}
		kh_value(hash, k) = (uint64_t)len<<32 | i;
		if (dret != '\n')
			while ((c = ks_getc(ks)) != '\n' && c != -1);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	fprintf(stderr, "[sam_header_read2] %d sequences loaded.\n", kh_size(hash));
	if (error) return 0;
	header = hash2header(hash);
	kh_destroy(ref, hash);
	return header;
}

static inline uint8_t *alloc_data(bam1_t *b, int size)
{
	if (b->m_data < size) {
		b->m_data = size;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	return b->data;
}

static inline void parse_error(int64_t n_lines, const char * __restrict msg)
{
	fprintf(stderr, "Parse error at line %lld: %s\n", (long long)n_lines, msg);
	abort();
}
static inline void append_text(bam_header_t *header, kstring_t *str)
{
	size_t x = header->l_text, y = header->l_text + str->l + 2; // 2 = 1 byte dret + 1 byte null
	kroundup32(x); kroundup32(y);
	if (x < y) 
    {
        header->n_text = y;
        header->text = (char*)realloc(header->text, y);
        if ( !header->text ) 
        {
            fprintf(stderr,"realloc failed to alloc %ld bytes\n", y);
            abort();
        }
    }
    // Sanity check
    if ( header->l_text+str->l+1 >= header->n_text )
    {
        fprintf(stderr,"append_text FIXME: %ld>=%ld, x=%ld,y=%ld\n",  header->l_text+str->l+1,(long)header->n_text,x,y);
        abort();
    }
	strncpy(header->text + header->l_text, str->s, str->l+1); // we cannot use strcpy() here.
	header->l_text += str->l + 1;
	header->text[header->l_text] = 0;
}

int sam_header_parse(bam_header_t *h)
{
	char **tmp;
	int i;
	free(h->target_len); free(h->target_name);
	h->n_targets = 0; h->target_len = 0; h->target_name = 0;
	if (h->l_text < 3) return 0;
	if (h->dict == 0) h->dict = sam_header_parse2(h->text);
	tmp = sam_header2list(h->dict, "SQ", "SN", &h->n_targets);
	if (h->n_targets == 0) return 0;
	h->target_name = calloc(h->n_targets, sizeof(void*));
	for (i = 0; i < h->n_targets; ++i)
		h->target_name[i] = strdup(tmp[i]);
	free(tmp);
	tmp = sam_header2list(h->dict, "SQ", "LN", &h->n_targets);
	h->target_len = calloc(h->n_targets, 4);
	for (i = 0; i < h->n_targets; ++i)
		h->target_len[i] = atoi(tmp[i]);
	free(tmp);
	return h->n_targets;
}

bam_header_t *sam_header_read(tamFile fp)
{
	int ret, dret;
	bam_header_t *header = bam_header_init();
	kstring_t *str = fp->str;
	while ((ret = ks_getuntil(fp->ks, KS_SEP_TAB, str, &dret)) >= 0 && str->s[0] == '@') { // skip header
		str->s[str->l] = dret; // note that str->s is NOT null terminated!!
		append_text(header, str);
		if (dret != '\n') {
			ret = ks_getuntil(fp->ks, '\n', str, &dret);
			str->s[str->l] = '\n'; // NOT null terminated!!
			append_text(header, str);
		}
		++fp->n_lines;
	}
	sam_header_parse(header);
	bam_init_header_hash(header);
	fp->is_first = 1;
	return header;
}

int alloc_sam(sam_msg *m, int xa_limit)	// for XA:Z: tags
{
	int i;
	if (m->sam_n >= xa_limit) return 1;
	if (m->sam_n >= m->sam_m) {
		m->sam_m <<= 1;
		if ((m->sam = (sam_t*)realloc(m->sam, m->sam_m * sizeof(sam_t))) == NULL) {fprintf(stderr, "[sam_parse] not enough mem.\n"); exit(-1);}
		for (i = m->sam_m/2; i < m->sam_m; ++i)
			m->sam[i].cigar_s = (kstring_t*)calloc(1, sizeof(kstring_t));;
	}
	return 0;
}

int sam_read1(tamFile fp, bam_header_t *header, sam_msg *m, int xa_limit)
{
	m->sam_n = 0;
	int ret, dret;
	kstring_t *str = fp->str;
	kstream_t *ks = fp->ks;

	if (fp->is_first) {
		fp->is_first = 0;
		ret = str->l;
	} else {
		do { // special consideration for empty lines
			ret = ks_getuntil(fp->ks, KS_SEP_TAB, str, &dret);
		} while (ret == 0);
	}
	if (ret < 0) return -1;
	++fp->n_lines;

	{ // name
        //printf("%s\t", str->s);
	}
	{ // flag
		long flag;
		char *s;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret);
		flag = strtol((char*)str->s, &s, 0);
		if (*s) { // not the end of the string
			flag = 0;
			for (s = str->s; *s; ++s)
				flag |= bam_char2flag_table[(int)*s];
		}
		if (flag & BAM_FUNMAP) { /*printf("%ld 0\n", flag);*/ ks_getuntil(ks, KS_SEP_LINE, str, &dret); return m->sam_n;}
		else if (flag & BAM_FREVERSE) {
			//printf("%ld -\t", flag);
			m->sam->nsrand = '-';
		}
		else {
			//printf("%ld +\t", flag);
			m->sam->nsrand = '+';
		}
		m->sam_n = 1;
	}
	{ // tid, pos, qual
		int tid, pos;//, qual;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret);
		tid = bam_get_tid(header, str->s);
        //printf("%d\t", tid);
		m->sam->chr = tid+1;
		if (tid < 0 && strcmp(str->s, "*")) {
			if (header->n_targets == 0) {
				fprintf(stderr, "[sam_read1] missing header? Abort!\n");
				exit(1);
			} else fprintf(stderr, "[sam_read1] reference '%s' is recognized as '*'.\n", str->s);
		}
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret);
		pos = isdigit(str->s[0])? atoi(str->s) : -1;    // 1-base
        //printf("%d\t", pos);
		m->sam->offset = pos;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret);
		//qual = isdigit(str->s[0])? atoi(str->s) : 0;
		if (ret < 0) return -2;
	}
	{ // cigar
		if (ks_getuntil(ks, KS_SEP_TAB, m->sam->cigar_s, &dret) < 0) return -3;
        //printf("%s\n", m->sam->cigar_s->s);
	}
	{ // mtid, mpos, isize
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret);
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret);
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret);
		if (ret < 0) return -4;
	}
	{ // seq and qual
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -5; // seq
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -6; // qual
	}
	{ // aux
		if (dret != '\n' && dret != '\r') { // aux exist
			char *s=0; int i, j;
			char *ts=0, *cigar=0, csrand;
			long long offset;
			int ed;
			while (ks_getuntil(ks, KS_SEP_TAB, str, &dret) >= 0) {
				if (!str->l) break;
				if (str->l < 6 || str->s[2] != ':' || str->s[4] != ':') parse_error(fp->n_lines, "missing colon in auxiliary data");
				if (str->s[0] == 'X' && str->s[1] == 'A' && str->s[3] == 'Z') 
				{
					//printf("\t\t%s\n", str->s);
					s = (char*)realloc(s, str->l - 4);
					ts = (char*)realloc(ts, str->l - 4);
					cigar = (char*)realloc(cigar, str->l - 4);
					j = 0;
					for (i = 5; i < str->l; ++i) {
						if (str->s[i] != ';')
							s[j++] = str->s[i];
						else {	//get a whole sam_msg
							if (alloc_sam(m, xa_limit)) { m->sam_n = 0; return 0; }
							s[j] = '\0';
							if ((sscanf(s, "%[^,],%c%lld,%[^,],%d", ts, &csrand, &offset, cigar, &ed)) != 5) return m->sam_n;

							m->sam[m->sam_n].chr = bam_get_tid(header, ts) + 1;
							m->sam[m->sam_n].nsrand = csrand;
							m->sam[m->sam_n].offset = offset;
							m->sam[m->sam_n].cigar_s->l = 0;
							kputs(cigar, m->sam[m->sam_n].cigar_s);
							++m->sam_n;
							j = 0;
						}
					}
					break;
				}
				if (dret == '\n' || dret == '\r') break;
			}
			if (s) free(s);
			if (ts) free(ts);
			if (cigar) free(cigar);
		}
	}
	return m->sam_n;
}

tamFile sam_open(const char *fn)
{
	tamFile fp;
	gzFile gzfp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "rb") : gzopen(fn, "rb");
	if (gzfp == 0) return 0;
	fp = (tamFile)calloc(1, sizeof(struct __tamFile_t));
	fp->str = (kstring_t*)calloc(1, sizeof(kstring_t));
	fp->fp = gzfp;
	fp->ks = ks_init(fp->fp);
	return fp;
}

void sam_close(tamFile fp)
{
	if (fp) {
		ks_destroy(fp->ks);
		gzclose(fp->fp);
		free(fp->str->s); free(fp->str);
		free(fp);
	}
}

// XXX end
