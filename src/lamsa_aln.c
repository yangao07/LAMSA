#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include "lamsa_aln.h"
#include "lamsa_dp_con.h"
#include "bwt.h"
#include "bntseq.h"
#include "frag_check.h"
#include "split_mapping.h"
#include "bwt_aln.h"
#include "kseq.h"
#include "kstring.h"
#include "gem_parse.h"

#define LINE_SIZE 65536
extern char lamsa_pg[1024];

int lamsa_aln_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   lamsa aln [options] <ref_prefix> <read.fa/fq>\n\n");

    fprintf(stderr, "Algorithm options:\n\n");

    fprintf(stderr, "    -t --thread    [INT]    Number of threads. [1]\n");
    //fprintf(stderr, "         -s [INT]      The seeding program, <gem(0)>, <bwa(1)> or <soap2-dp(2)>. [0]\n");
	fprintf(stderr, "    -l --seed-len  [INT]    Length of seeding fragments. [%d]\n", SEED_LEN);
	fprintf(stderr, "    -i --seed-inv  [INT]    Distance between neighboring seeding fragments. [%d]\n", SEED_INTERVAL);
	fprintf(stderr, "    -p --max-loci  [INT]    Maximum allowed number of seeding fragments' hits. [%d]\n", SEED_PER_LOCI);
    fprintf(stderr, "    -V --SV-len    [INT]    Expected maximum length of SV. [%d]\n", SV_MAX_LEN);
    fprintf(stderr, "    -v --ovlp-rat  [FLOAT]  Minimum overlapping ratio to cluster two skeletons or alignment records. (0~1) [%.1f]\n", OVLP_RAT);
    fprintf(stderr, "    -s --max-skel  [INT]    Maximum number of skeletons that are reserved in a cluster. [%d]\n", MAX_SKEL);
    fprintf(stderr, "    -R --max-reg   [INT]    Maximum allowed length of unaligned read part to trigger a bwt-based query. [%d]\n", MAX_BWT_REG);
    fprintf(stderr, "    -k --bwt-kmer  [INT]    Length of BWT-seed. [%d]\n\n", BWT_KMER);
	//fprintf(stderr, "         -d [INT]      The maximum number of seed's locations for first round's DP. [2]\n\n");
    

    fprintf(stderr, "Scoring options:\n\n");

    fprintf(stderr, "    -m --match-sc  [INT]    Match score for SW-alignment. [%d]\n", MAT_SCORE);
    fprintf(stderr, "    -M --mis-pen   [INT]    Mismatch penalty for SW-alignment. [%d]\n", MIS_PEN);
    fprintf(stderr, "    -O --open-pen  [INT]    Gap open penalty for SW-alignment. [%d]\n", OPEN_PEN);
    fprintf(stderr, "    -E --ext-pen   [INT]    Gap extension penalty for SW-alignment. [%d]\n\n", EXT_PEN);

    fprintf(stderr, "Output options:\n\n");

    fprintf(stderr, "    -r --max-out   [INT]    Maximum number of output records for a specific split read region. [%d]\n", RES_MAX_N);
    fprintf(stderr, "    -g --gap-split [INT]    Minimum length of gap that causes a split-alignment. [%d]\n", SPLIT_ALN_LEN);
    fprintf(stderr, "    -S --soft-clip          Use soft clipping for supplementary alignment. [false]\n");
    fprintf(stderr, "    -C --comment            Append FASTQ comment to SAM output. [false]\n");
    fprintf(stderr, "    -o --output    [STR]    Output file (SAM format). [stdout]\n\n");

    //fprintf(stderr, "         -N            Do NOT excute seeding program, when seeds' alignment result existed already.\n");
    //fprintf(stderr, "         -I            Seed information file has already existed.\n");
    //fprintf(stderr, "         -A [STR]      The seeds' alignment result. When '-N' is used. [\"seed_prefix.out.0\"]\n");
    
    fprintf(stderr, "    -h --help               Print this short usage.\n");
    fprintf(stderr, "    -H --HELP               Print a detailed usage.\n\n");
    return 1;
}

int lamsa_aln_de_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   lamsa aln [options] <ref_prefix> <read.fa/fq>\n\n");

    fprintf(stderr, "Algorithm options:\n\n");

    fprintf(stderr, "    -t --thread    [INT]    Number of threads. [1]\n");
    //fprintf(stderr, "         -s [INT]      The seeding program, <gem(0)>, <bwa(1)> or <soap2-dp(2)>. [0]\n");
	fprintf(stderr, "    -l --seed-len  [INT]    Length of seeding fragments. Moreover, LAMSA splits the read into a\n");
    fprintf(stderr, "                            series of -l bp long fragments, and employs NGS aligner to generate the\n");
    fprintf(stderr, "                            approximate matches of the fragments. [50]\n");
	fprintf(stderr, "    -i --seed-inv  [INT]    Distance between neighboring seeding fragments, LAMSA extracts seeding\n");
    fprintf(stderr, "                            fragments starting at every -i bp of the read. [%d]\n", SEED_INTERVAL);
	fprintf(stderr, "    -p --max-loci  [INT]    Maximum allowed number of hits. If a seeding fragment has more than -p\n");
    fprintf(stderr, "                            approximate matches, LAMSA would consider the seed is too repetitive, and\n");
    fprintf(stderr, "                            discard all the matches. [%d]\n", SEED_PER_LOCI);
    fprintf(stderr, "    -V --SV-len    [INT]    Expected maximum length of SV. If the genomic distance of two seeding\n");
    fprintf(stderr, "                            fragment is longer than -V bp, they cannot be connected to build a legal\n");
    fprintf(stderr, "                            edge in the sparse dynamic programming process. [%d]\n", SV_MAX_LEN);
    fprintf(stderr, "    -v --ovlp-rat  [FLOAT]  Minimum overlapping ratio to cluster two skeletons or alignment records.\n");
    fprintf(stderr, "                            (0~1) [%.1f]\n", OVLP_RAT);
    fprintf(stderr, "    -s --max-skel  [INT]    Maximum number of skeletons that are reserved in a cluster for a specific\n");
    fprintf(stderr, "                            read region. For a specific region of read, LAMSA reserves the top -s\n");
    fprintf(stderr, "                            skeletons. These skeletons are used to generate best and alternative\n");
    fprintf(stderr, "                            alignment records. [%d] \n", MAX_SKEL);
    fprintf(stderr, "    -R --max-reg   [INT]    Maximum allowed length of unaligned read part to trigger a bwt-based query.\n");
    fprintf(stderr, "                            For a read part being unaligned after all the sparse dynamic programming\n");
    fprintf(stderr, "                            process-based split alignment, if it is longer than -R bp, LAMSA will not\n");
    fprintf(stderr, "                            further process it; otherwise, LAMSA will query the exact matches of the\n");
    fprintf(stderr, "                            kmers of the unaligned part as hits to further align the read part. [%d]\n", MAX_BWT_REG);
    fprintf(stderr, "    -k --bwt-kmer  [INT]    Length of BWT-seed. For the unaligned read part shorter than -R bp, LAMSA\n");
    fprintf(stderr, "                            will extract all its -k bp tokens and query their exact match as hits. [%d]]\n\n", BWT_KMER);
	//fprintf(stderr, "         -d [INT]      The maximum number of seed's locations for first round's DP. [2]\n\n");
    
    fprintf(stderr, "Scoring options:\n\n");

    fprintf(stderr, "    -m --match-sc  [INT]    Match score for SW-alignment. [%d]\n", MAT_SCORE);
    fprintf(stderr, "    -M --mis-pen   [INT]    Mismatch penalty for SW-alignment. [%d]\n", MIS_PEN);
    fprintf(stderr, "    -O --open-pen  [INT]    Gap open penalty for SW-alignment. [%d]\n", OPEN_PEN);
    fprintf(stderr, "    -E --ext-pen   [INT]    Gap extension penalty for SW-alignment. A gap of length k costs O + k*E\n");
    fprintf(stderr, "                            (i.e. -O is for opening a zero-length gap). [%d]\n\n", EXT_PEN);

    fprintf(stderr, "Output options:\n\n");
    fprintf(stderr, "    -r --max-out   [INT]    Maximum number of output records for a specific split read region. For a\n");
    fprintf(stderr, "                            specific region, LAMSA reserves the top -r alignment records. The record\n");
    fprintf(stderr, "                            with highest alignment score is considered as best alignment, others are\n");
    fprintf(stderr, "                            considered as alternative alignments. If the score of an alternative\n");
    fprintf(stderr, "                            alignment is less than half of the best alignment, it will not be output.\n");
    fprintf(stderr, "                            [%d]\n", RES_MAX_N);
    fprintf(stderr, "    -g --gap-split [INT]    Minimum length of gap that causes a split-alignment. To avoid generating\n");
    fprintf(stderr, "                            insertion(I) or deletion(D) longer than -g bp in the SAM cigar. [%d]\n", SPLIT_ALN_LEN);
    fprintf(stderr, "    -S --soft-clip          Use soft clipping for supplementary alignment. It is strongly recommended\n");
    fprintf(stderr, "                            to turn off this option to reduce the redundancy of output when mapping\n");
    fprintf(stderr, "                            relatively long reads. [false]\n");
    fprintf(stderr, "    -C --comment            Append FASTQ comment to SAM output. [false]\n");
    fprintf(stderr, "    -o --output    [STR]    Output file (SAM format). [stdout]\n\n");

    //fprintf(stderr, "         -N            Do NOT excute seeding program, when seeds' alignment result existed already.\n");
    //fprintf(stderr, "         -I            Seed information file has already existed.\n");
    //fprintf(stderr, "         -A [STR]      The seeds' alignment result. When '-N' is used. [\"seed_prefix.out.0\"]\n");
    
    fprintf(stderr, "    -h --help               Print a short usage.\n");
    fprintf(stderr, "    -H --HELP               Print this detailed usage.\n\n");
    return 1;
}    

int f_BCC_score_table[10] = {
      1,  //F_MATCH
      1,  //F_SPLIT_MATCH
      1,  //F_MISMATCH
      1,  //F_LONG_MISMATCH
     -3,  //F_INSERT
     -3,  //F_DELETE
     -3,  //F_CHR_DIF
     -3,  //F_REVERSE
     -6,  //F_UNCONNECT
     -6,  //F_UNMATCH
};

// for debug
char READ_NAME[1024];

seed_msg *seed_init_msg(void)
{
    seed_msg *msg = (seed_msg*)malloc(sizeof(seed_msg));

    msg->read_all = 0;
    msg->read_m = READ_MAX_NUM;

    msg->seed_all = (int *)calloc(READ_MAX_NUM, sizeof(int));
    msg->read_len = (int *)malloc(READ_MAX_NUM * sizeof(int));
    /*
       msg->read_name = (char **)malloc(READ_MAX_NUM * sizeof(char*));
       int i;
       for (i = 0; i < READ_MAX_NUM; ++i)
        msg->read_name[i] = (char*)malloc(1024 * sizeof(char));
    */
    msg->last_len = (int *)malloc(READ_MAX_NUM * sizeof(int));
    msg->seed_max = 0;
    msg->read_max_len = 0;

    return msg;
}

void seed_free_msg(seed_msg *msg)
{
    free(msg->seed_all);
    free(msg->read_len);
    free(msg->last_len);
    free(msg);
}

int split_seed(const char *prefix, lamsa_aln_para AP, seed_msg *s_msg)
{
    gzFile infp;
    kseq_t *seq;
    char out_f[1024], seed_head[1024], seed_seq[1024], seed_info[1024];
    FILE *outfp, *infofp;
    int m_read, seed_all, i;
    void *new_p;

    if ((infp = gzopen(prefix, "r")) == NULL) {
        fprintf(stderr, "[lamsa_aln] Can't open read file %s\n", prefix); exit(1);
    }
    seq = kseq_init(infp);

    strcpy(out_f, prefix); strcat(out_f, ".seed");
    if ((outfp = fopen(out_f, "w")) == NULL) {
        fprintf(stderr, "[lamsa_aln] Can't open seed file %s\n", out_f); exit(1);
    }
    strcpy(seed_info, prefix); strcat(seed_info, ".seed.info");
    if ((infofp = fopen(seed_info, "w")) == NULL) {
        fprintf(stderr, "[lamsa_aln] Can't open seed info file %s\n", seed_info); exit(1);
    }

    fprintf(stderr, "[lamsa_aln] Spliting seed ... ");

    int seed_len=AP.seed_len, seed_inv=AP.seed_inv;
    m_read = s_msg->read_m;
    while (kseq_read(seq) >= 0)
    {
        seed_all = (1+ (seq->seq.l - seed_len) / (seed_len + seed_inv));
        seed_seq[seed_len] = '\n';
        if (seed_all > s_msg->seed_max) s_msg->seed_max = seed_all;
        if ((int)seq->seq.l > s_msg->read_max_len) s_msg->read_max_len = seq->seq.l;
        if (s_msg->read_all == m_read-1)
        {
            m_read <<= 1;
            if ((new_p = (int*)realloc(s_msg->seed_all, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_all);
                fprintf(stderr, "\n[lamsa_aln] Can't allocate more memory for seed_all[].\n"); exit(1);
            }
            s_msg->seed_all = (int*)new_p;
            if ((new_p = (int*)realloc(s_msg->last_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->last_len);
                fprintf(stderr, "\n[lamsa_aln] Can't allocate more memory for last_len[].\n"); exit(1);
            }
            s_msg->last_len = (int*)new_p;
            if ((new_p = (int*)realloc(s_msg->read_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_len);
                fprintf(stderr, "\n[lamsa_aln] Can't allocate more memory for read_len[].\n"); exit(1);
            }
            s_msg->read_len = (int*)new_p;
        }
        ++s_msg->read_all;
        s_msg->seed_all[s_msg->read_all] = s_msg->seed_all[s_msg->read_all-1] + seed_all;
        s_msg->last_len[s_msg->read_all] = seq->seq.l - seed_all * seed_len - (seed_all-1) * seed_inv;
        s_msg->read_len[s_msg->read_all] = seq->seq.l;
        s_msg->read_m = m_read;

        for (i = 0; i < seed_all; ++i)
        {
            sprintf(seed_head, ">%s_%d:%d\n", seq->name.s, i, i*(seed_len+seed_inv));
            strncpy(seed_seq, seq->seq.s+i*(seed_len+seed_inv), seed_len);
            seed_seq[seed_len+1] = '\0';
            fputs(seed_head, outfp);
            fputs(seed_seq, outfp);
        }
        fprintf(infofp, "%s %d %d %d\n", seq->name.s, seed_all, s_msg->last_len[s_msg->read_all], (int)seq->seq.l);
    }

    fprintf(stderr, "done!\n");
    gzclose(infp);
    fclose(outfp);
    fclose(infofp);
    kseq_destroy(seq);

    return 0;
}

int split_seed_info(const char *prefix, lamsa_aln_para AP, seed_msg *s_msg)
{
    char seed_info[1024];
    char read_name[1024];
    FILE *infofp;
    int m_read, seed_all, last_len, len, n;
    int seed_len=AP.seed_len, seed_inv=AP.seed_inv;
    void *new_p;

    strcpy(seed_info, prefix); strcat(seed_info, ".seed.info");
    if ((infofp = fopen(seed_info, "r")) == NULL)
    {
        fprintf(stderr, "[split seed] Can't open %s.\n", seed_info); exit(1);
    }
    m_read = s_msg->read_m;
    fprintf(stderr, "[lamsa_aln] Parsing seeds' information ... ");
    
    while ((n = fscanf(infofp, "%s %d %d %d", read_name, &seed_all, &last_len, &len)) != EOF)
    {
        if (n != 4)
        {
            fprintf(stderr, "\n[split seed] INFO file error.[2]\n"); exit(1);
        }
        if (seed_all > s_msg->seed_max) s_msg->seed_max = seed_all;
        if (len > s_msg->read_max_len) s_msg->read_max_len = len;
        if (s_msg->read_all == m_read-1)
        {
            m_read <<= 1;
            if ((new_p = (int*)realloc(s_msg->seed_all, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_all);
                fprintf(stderr, "\n[lamsa aln] Can't allocate more memory for seed_all[].\n"); exit(1);
            }
            s_msg->seed_all = (int*)new_p;
            if ((new_p = (int*)realloc(s_msg->last_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->last_len);
                fprintf(stderr, "\n[lamsa aln] Can't allocate more memory for last_len[].\n"); exit(1);
            }
            s_msg->last_len = (int*)new_p;
            if ((new_p = (int*)realloc(s_msg->read_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_len);
                fprintf(stderr, "\n[lamsa aln] Can't allocate more memory for read_len[].\n"); exit(1);
            }
            s_msg->read_len = (int*)new_p;
        }
        ++s_msg->read_all;
        s_msg->seed_all[s_msg->read_all] = s_msg->seed_all[s_msg->read_all-1] + seed_all;
        s_msg->last_len[s_msg->read_all] = last_len;
        s_msg->read_len[s_msg->read_all] = len;
        s_msg->read_m = m_read;
        if (last_len != len - seed_all * seed_len - (seed_all-1)*seed_inv)
        {
            fprintf(stderr, "\n%s %d %d %d %d %d", read_name, seed_all, last_len, len, seed_all, seed_inv);
            fprintf(stderr, "\n[split seed] INFO file error.[3]\n"); exit(1);
        }
    }

    fprintf(stderr, "done!\n");
    fclose(infofp);
    return 0;
}

aln_res *aln_init_res(int l_m, int n, int XA_max)
{
    int i, j, k;
    aln_res *a_res = (aln_res*)calloc(n, sizeof(aln_res));
    for (k = 0; k < n; ++k) {
        aln_res *p = a_res+k;
        p->l_m = l_m, p->l_n = 0;
        p->la = (line_aln_res*)malloc(l_m * sizeof(line_aln_res));
        //p->read_len = read_len;
        for (i = 0; i < l_m; ++i) {
            p->la[i].res_m = 10, p->la[i].cur_res_n = 0;
            p->la[i].res = (res_t*)malloc(10 * sizeof(res_t));
            for (j = 0; j < 10; ++j) {
                p->la[i].res[j].c_m = 100;
                p->la[i].res[j].cigar = (cigar32_t*)malloc(100 * sizeof(cigar32_t));
                p->la[i].res[j].cigar_len = 0;
            }
            p->la[i].tol_score = p->la[i].tol_NM = 0;
            // XA_res
            p->la[i].XA_m = XA_max;
            p->la[i].XA_n = 0;
            p->la[i].XA = (res_t**)malloc(XA_max * sizeof(res_t*));
        }
    }
    return a_res;
}

void aln_reset_res(aln_res *a_res, int n, int read_len)
{
    int i, j, k;
    for (k = 0; k < n; ++k) {
        aln_res *p = a_res+k;
        p->l_n = 0;
        p->read_len = read_len;
        for (i = 0; i < p->l_m; ++i) {
            p->la[i].cur_res_n = 0;
            for (j = 0; j < p->la[i].res_m; ++j) {
                p->la[i].res[j].cigar_len = 0;
            }
            p->la[i].tol_score = p->la[i].tol_NM = 0;
            // XA_res
            p->la[i].XA_n = 0;
        }
    }
}

void aln_reloc_res(aln_res *a_res, int line_n, int XA_m)
{
    int i, j;
    if ((a_res->la = (line_aln_res*)realloc(a_res->la, line_n * sizeof(line_aln_res))) == NULL) {
        fprintf(stderr, "[aln_reloc_res] Not enough memory.\n"); exit(0);
    }
    for (i = a_res->l_m; i < line_n; ++i) {
        a_res->la[i].res_m = 10, a_res->la[i].cur_res_n = 0;
        a_res->la[i].res = (res_t*)malloc(10 * sizeof(res_t));
        for (j = 0; j < 10; ++j) {
            a_res->la[i].res[j].c_m = 100;
            a_res->la[i].res[j].cigar = (cigar32_t*)malloc(100 * sizeof(cigar32_t));
            a_res->la[i].res[j].cigar_len = 0;
        }
        a_res->la[i].tol_score = a_res->la[i].tol_NM = 0;
        //XA
        a_res->la[i].XA_m = XA_m;
        a_res->la[i].XA_n = 0;
        a_res->la[i].XA = (res_t**)malloc(XA_m * sizeof(res_t*));
    }
    a_res->l_m = line_n;
}

void aln_res_free(aln_res *res, int n)
{
    int i,j,k;
    for (k = 0; k < n; ++k) {
        aln_res *p = res+k;
        for (i = 0; i < p->l_m; ++i) {
            for (j = 0; j < p->la[i].res_m; ++j) {
                free(p->la[i].res[j].cigar);
            }
            free(p->la[i].res); free(p->la[i].XA);
        }
        free(p->la);
    }
    free(res);
}

aln_reg *aln_init_reg(int read_len)
{
    aln_reg *reg = (aln_reg*)malloc(sizeof(aln_reg));
    reg->reg_n = 0; reg->reg_m = 1;
    reg->read_len = read_len;
    reg->reg = (reg_t*)malloc(sizeof(reg_t));
    reg->reg->beg_n = reg->reg->end_n = 0;
    reg->reg->beg_m = reg->reg->end_m = 10; 
    reg->reg->ref_beg = (reg_b*)malloc(10 * sizeof(reg_b));
    reg->reg->ref_end = (reg_b*)malloc(10 * sizeof(reg_b));

    return reg;
}

void aln_free_reg(aln_reg *reg) { 
    int i;
    for (i = 0; i < reg->reg_m; ++i) {
        free(reg->reg[i].ref_beg); free(reg->reg[i].ref_end);
    }
    free(reg->reg); free(reg); 
}

int reg_comp(const void *a, const void *b) { return (*(reg_t*)a).beg - (*(reg_t*)b).beg; }
void aln_sort_reg(aln_reg *a_reg) { qsort(a_reg->reg, a_reg->reg_n, sizeof(reg_t), reg_comp); }

void push_reg_b(reg_t *reg1, int beg_n, reg_b beg[], int end_n, reg_b end[]) {
    int i;
    for (i = 0; i < beg_n; ++i) {
        if (reg1->beg_n == reg1->beg_m) {
            reg1->beg_m <<= 1;
            reg1->ref_beg = (reg_b*)realloc(reg1->ref_beg, reg1->beg_m * sizeof(reg_b));
        }
        reg1->ref_beg[reg1->beg_n] = beg[i];
        reg1->beg_n++;
    }
    for (i = 0; i < end_n; ++i) {
        if (reg1->end_n == reg1->end_m) {
            reg1->end_m <<= 1;
            reg1->ref_end = (reg_b*)realloc(reg1->ref_end, reg1->end_m * sizeof(reg_b));
        }
        reg1->ref_end[reg1->end_n] = end[i];
        reg1->end_n++;
    }
}

void aln_merg_reg(aln_reg *a_reg, int thd) {
    //if (a_reg->reg_n == 0)return;
    int cur_i = 0, i;
    for (i = 1; i < a_reg->reg_n; ++i) {
		// merge when the distance between two regs is smaller than the 'thd'
        if(a_reg->reg[i].beg - a_reg->reg[cur_i].end - 1 < thd) {
		    if(a_reg->reg[i].end > a_reg->reg[cur_i].end) a_reg->reg[cur_i].end = a_reg->reg[i].end;
            push_reg_b(a_reg->reg+cur_i, a_reg->reg[i].beg_n, a_reg->reg[i].ref_beg,a_reg->reg[i].end_n, a_reg->reg[i].ref_end);
		} else {
            cur_i++;
            if (cur_i != i) {
                a_reg->reg[cur_i].beg = a_reg->reg[i].beg; 
                a_reg->reg[cur_i].end = a_reg->reg[i].end;
                a_reg->reg[cur_i].beg_n = a_reg->reg[cur_i].end_n = 0;
                push_reg_b(a_reg->reg+cur_i, a_reg->reg[i].beg_n, a_reg->reg[i].ref_beg,a_reg->reg[i].end_n, a_reg->reg[i].ref_end);
            }
        }
    }
    a_reg->reg_n = cur_i+1;
}

void reg_reloc(reg_t **r, int r_n, int r_m) {
    *r = (reg_t*)realloc(*r, r_m * sizeof(reg_t));
    int i;
    for (i = r_n; i < r_m; ++i) {
        (*r)[i].end_n = (*r)[i].beg_n = 0;
        (*r)[i].end_m = (*r)[i].beg_m = 10;
        (*r)[i].ref_beg = (reg_b*)malloc(10 * sizeof(reg_b));
        (*r)[i].ref_end = (reg_b*)malloc(10 * sizeof(reg_b));
    }
}

void push_reg(aln_reg *reg, int beg, int end, int beg_n, reg_b ref_beg[], int end_n, reg_b ref_end[])
{
    int r_n = reg->reg_n;
    if (r_n == reg->reg_m) {
        reg->reg_m <<= 1;
        reg_reloc(&(reg->reg), r_n, reg->reg_m);
    }
    reg->reg[r_n].beg = beg;
    reg->reg[r_n].end = end;
    reg->reg[r_n].beg_n = reg->reg[r_n].end_n = 0;
    push_reg_b(reg->reg+r_n, beg_n, ref_beg, end_n, ref_end);
    reg->reg_n++;
}
/*typedef struct {
	int refid, is_rev;
	uint64_t ref_beg, ref_end;
    int beg, end;
} reg_t;*/
// dump remain-regions that are longger than 'thd'
int get_remain_reg(aln_reg *a_reg, aln_reg *remain_reg, lamsa_aln_para AP, int reg_thd)
{
    if (a_reg->reg_n == 0) {
        if (AP.bwt_seed_len < a_reg->read_len && a_reg->read_len <= reg_thd) {
            push_reg(remain_reg, 1, a_reg->read_len, 0, 0, 0, 0);
            return 1;
        } else return 0;
    }
    aln_sort_reg(a_reg); aln_merg_reg(a_reg, AP.bwt_seed_len);
    int i; 
    if (a_reg->reg[0].beg > AP.bwt_seed_len && a_reg->reg[0].beg-1 <= reg_thd)
        push_reg(remain_reg, 1, a_reg->reg[0].beg-1, 0, 0, a_reg->reg[0].beg_n, a_reg->reg[0].ref_beg);
    for (i = 1; i < a_reg->reg_n; ++i) {
        if (a_reg->reg[i].beg - a_reg->reg[i-1].end > AP.bwt_seed_len && a_reg->reg[i].beg-1-a_reg->reg[i-1].end <= reg_thd)
            push_reg(remain_reg, a_reg->reg[i-1].end+1, a_reg->reg[i].beg-1, a_reg->reg[i-1].end_n, a_reg->reg[i-1].ref_end, a_reg->reg[i].beg_n, a_reg->reg[i].ref_beg);
    }
	if (a_reg->read_len - a_reg->reg[i-1].end > AP.bwt_seed_len && a_reg->read_len-a_reg->reg[i-1].end <= reg_thd) 
        push_reg(remain_reg, a_reg->reg[i-1].end+1, a_reg->read_len, a_reg->reg[i-1].end_n, a_reg->reg[i-1].ref_end, 0, 0);
    return remain_reg->reg_n;
}

void push_reg_res(aln_reg *reg, res_t *res)
{
    int r_n = reg->reg_n;
    if (r_n == reg->reg_m) {
        reg->reg_m <<= 1;
        reg_reloc(&(reg->reg), r_n, reg->reg_m);
    }
    reg->reg[r_n].ref_beg[0].chr = reg->reg[r_n].ref_end[0].chr = res->chr;
	reg->reg[r_n].ref_beg[0].is_rev = reg->reg[r_n].ref_end[0].is_rev = 1-res->nsrand;
    if (res->nsrand == 1) { // '+'
        reg->reg[r_n].beg = (res->cigar[0] & 0xf) == CSOFT_CLIP ? ((res->cigar[0]>>4)+1):1;
        reg->reg[r_n].end = (res->cigar[res->cigar_len-1] & 0xf) == CSOFT_CLIP ? (reg->read_len - (res->cigar[res->cigar_len-1]>>4)):reg->read_len;
		reg->reg[r_n].ref_beg[0].ref_pos = res->offset;
		reg->reg[r_n].ref_end[0].ref_pos = res->offset+refInCigar(res->cigar, res->cigar_len)-1;
    } else { // '-'
        reg->reg[r_n].beg = (res->cigar[res->cigar_len-1] & 0xf) == CSOFT_CLIP ? ((res->cigar[res->cigar_len-1]>>4) + 1):1;
        reg->reg[r_n].end = (res->cigar[0] & 0xf) == CSOFT_CLIP ? (reg->read_len - (res->cigar[0]>>4)):reg->read_len;
		reg->reg[r_n].ref_end[0].ref_pos = res->offset;
		reg->reg[r_n].ref_beg[0].ref_pos = res->offset+refInCigar(res->cigar, res->cigar_len)-1;
    }
    res->reg_beg = reg->reg[r_n].beg;
    res->reg_end = reg->reg[r_n].end;
    reg->reg[r_n].beg_n = reg->reg[r_n].end_n = 1;
    reg->reg_n++;
}

void get_reg(aln_res *res, aln_reg *reg)
{
    int i, j;
    for (i = 0; i < res->l_n; ++i) {
        if (res->la[i].tol_score < 0) continue;
        for (j = 0; j <= res->la[i].cur_res_n; ++j)
            push_reg_res(reg, res->la[i].res+j);
    }
}

int get_cover_res(aln_reg *reg, aln_res *res, int qua_i, int *cov_qua_i, qua_node *qua, int head[], int head_n, float ovlp_r)
{
    extern float cover_rate(int s1, int e1, int s2, int e2);
    int reg_i, i, j, res_i, l_i;
    int _res_i = qua[qua_i].x, _l_i = qua[qua_i].y, _r_i; // new res

    for (_r_i = 0; _r_i <= (res+_res_i)->la[_l_i].cur_res_n; ++_r_i) {
        res_t *_r = (res+_res_i)->la[_l_i].res+_r_i;
        reg_i = 0;
        for (i = 0; i < head_n; ++i) { // existing ress and regs
            res_i = qua[head[i]].x;
            l_i = qua[head[i]].y;
            for (j = 0; j <= (res+res_i)->la[l_i].cur_res_n; ++j) {
                if (cover_rate(reg->reg[reg_i].beg, reg->reg[reg_i].end, _r->reg_beg, _r->reg_end) >= ovlp_r) {
                    *cov_qua_i = head[i];
                    return 1; // cover
                }
                reg_i++;
            }
        }
    }
    return 0;// Not cover
}


int res_comp(const void*a,const void*b) 
{ 
    if ((*(qua_node*)b).a-(*(qua_node*)a).a)
        return (*(qua_node*)b).a-(*(qua_node*)a).a; 
    else return (*(qua_node*)b).b-(*(qua_node*)a).b; 
}

float get_cov_f(aln_res *res, aln_reg *reg)
{
    float cov_f; int i, cov=0;

    reg->reg_n = 0;
    for (i = 0; i < 3; ++i) get_reg(res+i, reg);
    aln_sort_reg(reg); aln_merg_reg(reg, 0);

    for (i = 0; i < reg->reg_n; ++i) 
        cov += (reg->reg[i].end-reg->reg[i].beg+1); 
    cov_f = (cov+0.0)/reg->read_len;
    return cov_f;
}

// rearrange the aln_res (merge and filter by the score and cover-region)
void rearr_aln_res(aln_res *res, int n, float ovlp_r)
{
    extern void copy_res(res_t *f, res_t *t);
    int a_i, i, j; aln_res *p;

    // sort by score;
    qua_node *qua = (qua_node*)malloc(10 * sizeof(qua_node));
    int qua_n = 0, qua_m = 10;
    int head[1024], head_n=0;
    
    for (a_i = 0; a_i < n; ++a_i) {
        p = res+a_i;
        for (i = 0; i < p->l_n; ++i) {
            if (p->la[i].tol_score < 0) { p->la[i].merg_msg = (line_node){0, -1};  continue; }
            if (qua_n == qua_m) {
                qua = (qua_node*)realloc(qua, qua_m * 2 * sizeof(qua_node));
                qua_m <<= 1;
            }
            qua[qua_n++] = (qua_node){a_i, i, p->la[i].tol_score, p->la[i].line_score};
        }
    }
    if (qua_n == 0) {free(qua); return;}

    qsort(qua, qua_n, sizeof(qua_node), res_comp);
    aln_reg *reg = aln_init_reg(res->read_len);
    // set the res with the biggest score: qua[0]
    for (i = 0; i <= (res+qua[0].x)->la[qua[0].y].cur_res_n; ++i) push_reg_res(reg, (res+qua[0].x)->la[qua[0].y].res+i);
    head[head_n++] = 0;
    // NOT merg
    (res+qua[0].x)->la[qua[0].y].merg_msg = (line_node){1, 0}; // NOT merge and ONLY best
    for (i = 0; i < qua_n; ++i) (res+qua[i].x)->la[qua[i].y].mapQ = 255; 

    int mapq_max = MAPQ_MAX * res->cov_f;
    
    int cov_qua_i;
    for (i = 1; i < qua_n; ++i) {
        if (!get_cover_res(reg, res, i, &cov_qua_i, qua, head, head_n, ovlp_r)) { // NOT cover
            // NOT merg
            for (j = 0; j <= (res+qua[i].x)->la[qua[i].y].cur_res_n; ++j) push_reg_res(reg, (res+qua[i].x)->la[qua[i].y].res+j);
            head[head_n++] = i;
            (res+qua[i].x)->la[qua[i].y].merg_msg = (line_node){1, 0}; // NOT merge and ONLY best
        } else if (qua[i].a > qua[cov_qua_i].a/2 
               && (res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].XA_n + (res+qua[i].x)->la[qua[i].y].cur_res_n 
               < (res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].XA_m) {
            // head_i.merg_flag = MERG_HEAD
            for (j = 0; j <= (res+qua[i].x)->la[qua[i].y].cur_res_n; ++j) {
                (res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].XA[(res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].XA_n] = (res+qua[i].x)->la[qua[i].y].res+j;
                (res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].XA_n++;
            }
            (res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].merg_msg.y = 1; // {1,1}, Merge, has alternative
            // set mapQ for primary
            uint8_t tmpQ = mapq_max * (qua[cov_qua_i].a - qua[i].a)/qua[cov_qua_i].a;
            if (tmpQ < ((res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].mapQ))
                (res+qua[cov_qua_i].x)->la[qua[cov_qua_i].y].mapQ = tmpQ;

            (res+qua[i].x)->la[qua[i].y].merg_msg = (line_node){2, 0}; // Merge, body 
        } else {
            // DUMP
            (res+qua[i].x)->la[qua[i].y].merg_msg = (line_node){0, -1}; // DUMP
        }
    }

    // re-set mapQ, based on head_n
    for (i = 0; i < qua_n; ++i) {
        if ((res+qua[i].x)->la[qua[i].y].merg_msg.x != 1) continue;
        if ((res+qua[i].x)->la[qua[i].y].mapQ == 255) {
            (res+qua[i].x)->la[qua[i].y].mapQ = mapq_max / head_n;
        } else (res+qua[i].x)->la[qua[i].y].mapQ /= head_n;
    }
    free(qua); aln_free_reg(reg);
}

aln_msg *aln_init_msg(int seed_max, int per_aln_m)
{
    aln_msg *msg;
    int i,j;
    msg = (aln_msg*)malloc(seed_max * sizeof(aln_msg));
    for (i = 0; i < seed_max; ++i)		//drop away seed whose number of alignments > PER_ALN_N
    {
        msg[i].read_id = -1;    // -> seed_id ?
        msg[i].n_aln = 0;
        msg[i].skip = 0;
        msg[i].at = (aln_t*)malloc(per_aln_m * sizeof(aln_t));  // per_aln_m
        for (j = 0; j < per_aln_m; ++j)
        {
            msg[i].at[j].cigar = (cigar32_t*)malloc(7 * sizeof(cigar32_t));//default value for 3-ed
            msg[i].at[j].cigar_len = 0;
            msg[i].at[j].cmax = 7;
            msg[i].at[j].bmax = 0;
        }
    }
    return msg;
}

void aln_free_msg(aln_msg *a_msg, int seed_max, int per_aln_m)	//a_msg[seed_max]
{
    int i,j;
    for (i = 0; i < seed_max; ++i)
    {
        for (j = 0; j < per_aln_m; ++j)
        {
            free(a_msg[i].at[j].cigar);
        }
        free(a_msg[i].at);
    }
    free(a_msg);
}

frag_dp_node ***fnode_alloc(int seed_m, int per_aln_m)
{
    int i, j;
    frag_dp_node ***f_node;
    f_node = (frag_dp_node***)malloc(sizeof(frag_dp_node**));
    (*f_node) = (frag_dp_node**)malloc(seed_m * sizeof(frag_dp_node*));
    for (i = 0; i < seed_m; ++i) {
        (*f_node)[i] = (frag_dp_node*)malloc(per_aln_m * sizeof(frag_dp_node));
        for (j = 0; j < per_aln_m; ++j) {
            (*f_node)[i][j].seed_i = i, (*f_node)[i][j].aln_i = j;
            (*f_node)[i][j].son_max = 100; // son_max 
            (*f_node)[i][j].son = (line_node*)calloc(100, sizeof(line_node));
        }
    }
    return f_node;
}

void fnode_init(frag_dp_node **f_node, int seed_m, int per_aln_m)
{
    int i, j;
    for (i = 0; i < seed_m; ++i) {
        for (j = 0; j < per_aln_m; ++j) {
            f_node[i][j].in_de = 0;
            f_node[i][j].son_n = 0;
            f_node[i][j].max_score = 0;
            f_node[i][j].trg_n = 0;
        }
    }
}

void fnode_free(frag_dp_node ***f_node, int seed_m, int per_aln_m)
{
    int i, j;
    for (i = 0; i < seed_m; ++i) {
        for (j = 0; j < per_aln_m; ++j)
            free((*f_node)[i][j].son);
        free((*f_node)[i]);
    }
    free(*f_node);
    free(f_node);
}

//MIDNSHP=XB
//0123456789
void setCigar(aln_msg *a_msg, int seed_i, int aln_i, char *s_cigar)
{
    int op;
    long x, bi, bd;
    char *s, *t;

    a_msg[seed_i].at[aln_i].cigar_len=0;
    bi = bd = 0;
    for (s = s_cigar; *s; ) {
        x = strtol(s, &t, 10);	
        if (x == 0) {
            fprintf(stderr, "\n[lamsa_aln] Cigar ERROR 1.\n"); fprintf(stderr, "%s\n",s); exit(1);
        }
        op = toupper(*t);
        switch (op)
        {
            case 'M':	op = CMATCH;	break;
            case 'I':	op = CINS;		bi += x;	break;
            case 'D':	op = CDEL;		bd += x;	break;
                        //case 'S':	op = CSOFT_CLIP;		bi += x;	break;
                        //case 'H':	op = CHARD_CLIP;		bi += x;	break;
            default:	fprintf(stderr, "\n[lamsa_aln] Cigar ERROR 2.\n"); exit(1); break;
        }
        if (a_msg[seed_i].at[aln_i].cigar_len == a_msg[seed_i].at[aln_i].cmax)
        {
            a_msg[seed_i].at[aln_i].cmax <<= 2 ;
            a_msg[seed_i].at[aln_i].cigar = (cigar32_t*)realloc(a_msg[seed_i].at[aln_i].cigar, a_msg[seed_i].at[aln_i].cmax * sizeof(cigar32_t));
        }
        a_msg[seed_i].at[aln_i].cigar[a_msg[seed_i].at[aln_i].cigar_len] = CIGAR_GEN(x, op);
        //modify variable directly OR use a auxiliary-variable
        ++a_msg[seed_i].at[aln_i].cigar_len;
        s = t+1;
    }
    a_msg[seed_i].at[aln_i].len_dif = (int)(bd - bi);
    a_msg[seed_i].at[aln_i].bmax = (int)(bd > bi ? bd : bi);
}

void set_cigar(aln_t *at, cigar_t *cigar)
{
    int i, bd=0, bi=0;
    at->cigar_len = 0;
    for (i = 0; i < cigar->cigar_n; ++i) {
        _push_cigar1(&(at->cigar), &(at->cigar_len), &(at->cmax), cigar->cigar[i]);
        if (((cigar->cigar[i]) & 0xf) == CINS) bi += (cigar->cigar[i] >> 4);
        else if (((cigar->cigar[i]) & 0xf) == CDEL) bd += (cigar->cigar[i] >> 4);
    }
    at->len_dif = bd - bi;
    at->bmax = bd > bi ? bd : bi;
}

void set_aln_msg(aln_msg *a_msg, int32_t read_x, int aln_y, int read_id, map_t map, bntseq_t *bns)
{
    a_msg[read_x-1].read_id = read_id;	
    a_msg[read_x-1].at[aln_y-1].chr = bns_get_rid(bns, map.chr);
    a_msg[read_x-1].at[aln_y-1].offset = map.offset;	//1-base
    a_msg[read_x-1].at[aln_y-1].nsrand = ((map.strand=='+')?1:-1);
    a_msg[read_x-1].at[aln_y-1].NM = map.NM;
    a_msg[read_x-1].n_aln = aln_y;
    set_cigar(a_msg[read_x-1].at+aln_y-1, map.cigar);
}

void init_aln_per_para(lamsa_aln_per_para *APP, seed_msg *s_msg, int read_n)
{
    APP->last_len = s_msg->last_len[read_n];
    APP->seed_all = s_msg->seed_all[read_n] - s_msg->seed_all[read_n-1];
}

typedef struct {
    //char *name;          // read name
    //uint8_t *seq;         // 0/1/2/3/4:A/C/G/T/N 
    //uint8_t *rseq;
    //char *qual;
    //int len;             // read length
    
    lamsa_aln_per_para *APP; // parameter for this read

    map_msg *m_msg;      // seeds' mapping infomation
                         // size: seed_all * sizeof(map_msg)

    aln_res *a_res;      // [3] alignment results
} lamsa_seq_t;            // whole_seqs_size = n_seqs * sizeof(lamsa_seq_t)

typedef struct {
    int tid;             // local id of thread
    lamsa_aln_para *AP;   // common parameter for all reads
    bwt_t *bwt;          // bwt index
    uint8_t *pac;        // pac of reference
    bntseq_t *bns;       // reference

    int n_seqs;
    lamsa_seq_t *lamsa_seqs;    // auxiliary: rseq, APP, m_msg, a_res
    kseq_t *w_seqs;       // whole seqs to be processed

    aln_msg *a_msg;      // seeds' mapping information
    frag_dp_node ***f_node; // DP nodes
    line_node **line, **_line;
    int *line_end, *_line_end;

    int line_n_max; int line_m;
    frag_msg **f_msg;     // alignment information of fragments of read

    uint32_t *hash_num;  // hash index
    uint64_t **hash_node;
} thread_aux_t;          // whole_aux_size = n_thread * sizeof(thread_aux_t)


int COUNT=0;
int THREAD_READ_I;
pthread_rwlock_t RWLOCK;
int lamsa_main_aln(thread_aux_t *aux)
{
    lamsa_aln_para *AP = aux->AP; bwt_t *bwt = aux->bwt; uint8_t *pac = aux->pac; bntseq_t *bns = aux->bns;
    aln_msg *a_msg = aux->a_msg; frag_dp_node ***f_node = aux->f_node;
    line_node **line = aux->line, **_line = aux->_line;
    int *line_end = aux->line_end; int *_line_end = aux->_line_end;
    frag_msg **f_msg = aux->f_msg;
    uint32_t *hash_num = aux->hash_num; uint64_t **hash_node = aux->hash_node;
    int line_n_max = aux->line_n_max, line_m = aux->line_m;

    int i = 0, j, k;

    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_READ_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= aux->n_seqs) break;
        lamsa_seq_t *la_seqs = aux->lamsa_seqs + i; lamsa_aln_per_para *APP = la_seqs->APP;
        kseq_t *seqs = aux->w_seqs+i;
        strcpy(READ_NAME, seqs->name.s);
        // set aln_msg
        int seed_out = 0;
        for (j = 0; j < APP->seed_all; ++j) {
            gem_map_msg(la_seqs->m_msg+j, AP->per_aln_m);
            if ((la_seqs->m_msg+j)->map_n > 0) seed_out++;
            APP->seed_out = seed_out;
            for (k = 0; k < (la_seqs->m_msg+j)->map_n; ++k)
                set_aln_msg(a_msg, seed_out, k+1, j+1, (la_seqs->m_msg+j)->map[k], bns);
        }
        // aln_res
        aln_reset_res(la_seqs->a_res, 3, seqs->seq.l);
        // aln_reg
        aln_reg *a_reg = aln_init_reg(seqs->seq.l);
        int line_n = frag_line_BCC(a_msg, f_msg, *APP, *AP, seqs, line, line_end, &line_m, f_node, _line, line_n_max);
        
        // char -> 0123
        uint8_t *bseq = (uint8_t*)malloc(seqs->seq.l * sizeof(uint8_t));
        for (j = 0; j < seqs->seq.l; ++j) bseq[j] = nst_nt4_table[(int)(seqs->seq.s[j])];
        uint8_t *rbseq=NULL;
        if (line_n > 0) {
            frag_check(a_msg, f_msg, la_seqs->a_res, bns, pac, bseq, &rbseq, *APP, *AP, seqs, line_n, &hash_num, &hash_node);
            get_reg(la_seqs->a_res, a_reg);
        }
        // remain region
        line_n = frag_line_remain(a_reg, a_msg, f_msg, *APP, *AP, seqs, line, line_end, &line_m, f_node, _line, _line_end, line_n_max);
        if (line_n > 0) {
            frag_check(a_msg, f_msg, la_seqs->a_res+1, bns, pac, bseq, &rbseq, *APP, *AP, seqs, line_n, &hash_num, &hash_node);
            get_reg(la_seqs->a_res+1, a_reg);
        }
        // bwt aln
        bwt_aln_remain(a_reg, la_seqs->a_res+2, bwt, bns, pac, bseq, &rbseq, *AP, seqs);
        get_reg(la_seqs->a_res+2, a_reg);

        la_seqs->a_res->cov_f = get_cov_f(la_seqs->a_res, a_reg);
        // rearrange a_res(0,1,2)
        rearr_aln_res(la_seqs->a_res, 3, AP->ovlp_rat);

        aln_free_reg(a_reg); if (rbseq) free(rbseq); free(bseq);

#ifdef __DEBUG__
        COUNT++; fprintf(stderr, "%16d reads have been aligned.\n", COUNT);
#else
        COUNT++;
        if (COUNT % 100000 == 0) fprintf(stderr, "%16d reads have been aligned.\n", COUNT);
#endif
    }
    aux->line_m = line_m;
    return 0;
}

static void *lamsa_thread_aln(void *aux)
{
    thread_aux_t *a = (thread_aux_t*)aux;
    lamsa_main_aln(a);
    return 0;
}

void bseq_reco(char *seq, int len)
{
    int i;
    for (i = 0; i < len>>1; ++i) {
        char tmp = seq[len-1-i];
        if (tmp < 4) tmp = 3 - tmp;
        seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3-seq[i];
        seq[i] = tmp;
    }
    if (len & 1) seq[i] = (seq[i] >= 4)? seq[i] : 3-seq[i];
}

lamsa_seq_t *lamsa_seq_init(int chunk_read_n, int seed_m, int x)
{
    int i;
    lamsa_seq_t *seqs = (lamsa_seq_t*)calloc(chunk_read_n, sizeof(lamsa_seq_t));
    if (seqs == NULL) { fprintf(stderr, "[lamsa_seq_init] Not enough memory.\n"); exit(1); }
    for (i = 0; i < chunk_read_n; ++i) {
        seqs[i].a_res = aln_init_res(1, 3, x); // aln_res, remain_res, bwt_remain_res
        seqs[i].APP = (lamsa_aln_per_para*)malloc(sizeof(lamsa_aln_per_para));
        seqs[i].m_msg = map_init_msg(seed_m);
    }
    return seqs;
}

// => seqs(name, seq, (map_msg,APP) * seed_num), n_seqs
// => seqs(aln_res) in thread
int lamsa_read_seq(lamsa_seq_t *la_seqs, kseq_t *read_seq_t, FILE *seed_mapfp, char *gem_line, int line_size, seed_msg *s_msg, int chunk_read_n)
{
    kseq_t *s = read_seq_t;
    lamsa_seq_t *p=NULL; int n = 0;
    while (kseq_read(s+n) >= 0)
    {
        p = &la_seqs[n++];

        // APP, seed_mapfp
        ++(s_msg->read_count);
        init_aln_per_para(p->APP, s_msg, s_msg->read_count);

        int seed_all = p->APP->seed_all;
        int seed_n = 0;
        while (seed_n < seed_all) {
            //if ((map_n = gem_map_read(seed_mapfp, p->m_msg+seed_n, AP.per_aln_m)) < 0) { fprintf(stderr, "[lamsa_read_seq] Seeds' GEM map-result do NOT match.\n"); exit(1); }
            if (gem_map_read(seed_mapfp, p->m_msg+seed_n, gem_line, line_size) < 0) { fprintf(stderr, "[lamsa_read_seq] Seeds' GEM map-result do NOT match.\n"); exit(1); }
            seed_n++;
        }
        if (n >= chunk_read_n) break;
    }
    return n;
}

void lamsa_free_read_seq(lamsa_seq_t *seqs, int n_seqs)
{
    int i;
    for (i = 0; i < n_seqs; ++i) {
        lamsa_seq_t *p = seqs+i;
        map_free_msg(p->m_msg, p->APP->seed_all);
        free(p->APP); 
        aln_res_free(p->a_res, 3);
    }
    free(seqs);
}

void aux_dp_init(thread_aux_t *aux, seed_msg *s_msg, lamsa_aln_para AP)
{
    aux->a_msg = aln_init_msg(s_msg->seed_max, AP.per_aln_m);  
    aux->f_node = fnode_alloc(s_msg->seed_max+2, AP.per_aln_m);

    int i, line_n_max = s_msg->seed_max * AP.per_aln_m, line_len_max = s_msg->seed_max;
    aux->line = (line_node**)malloc(line_n_max * sizeof(line_node*));
    aux->_line = (line_node**)malloc(line_n_max * sizeof(line_node*));
    aux->line_end = (int*)malloc(line_n_max * sizeof(int));
    aux->_line_end = (int*)malloc(line_n_max * sizeof(int));
    for (i = 0; i < line_n_max; ++i) {
        aux->line[i] = (line_node*)malloc((line_len_max+L_EXTRA) * sizeof(line_node));
        aux->_line[i] = (line_node*)malloc((line_len_max+L_EXTRA) * sizeof(line_node));
    }
    aux->line_n_max = line_n_max;
    aux->line_m = 1;

    aux->f_msg = (frag_msg**)malloc(sizeof(frag_msg*));
    *(aux->f_msg) = (frag_msg*)malloc(sizeof(frag_msg));
    frag_init_msg(*(aux->f_msg), s_msg->seed_max);

    aux->hash_num = (uint32_t*)calloc(pow(NT_N, AP.hash_key_len), sizeof(uint32_t));
    aux->hash_node = (uint64_t**)calloc(pow(NT_N, AP.hash_key_len), sizeof(uint64_t*));
}

void aux_dp_free(thread_aux_t *aux, seed_msg *s_msg, lamsa_aln_para *AP)
{
    aln_free_msg(aux->a_msg, s_msg->seed_max, AP->per_aln_m);
    fnode_free(aux->f_node, s_msg->seed_max+2, AP->per_aln_m);
    int i, line_n_max = s_msg->seed_max * AP->per_aln_m;
    for (i = 0; i < line_n_max; ++i) {
        free(aux->line[i]); free(aux->_line[i]);
    } free(aux->line); free(aux->_line); free(aux->line_end); free(aux->_line_end);
    
    for (i = 0; i < pow(NT_N, AP->hash_key_len); ++i) free(aux->hash_node[i]);
    free(aux->hash_node); free(aux->hash_num);

    frag_free_msg(*(aux->f_msg), aux->line_m); free(aux->f_msg);
}

//merg_msg: {1, 0} -> NOT merged or ONLY best
//          {1, n} -> merged, best, has n alternative res
//          {2, i} -> merged, alternative and best is `i`
//          {0,-1} -> dumped
void aln_res_output(lamsa_aln_para AP, aln_res *res, int res_n, char *name, char *seq, char *qual, bntseq_t *bns)
{
    FILE *outp = AP.outp;
    int i, j, n, l, all=0;
    kstring_t sam_str; memset(&sam_str, 0, sizeof(kstring_t));

    int prim_flag = 0, sam_flag;
    for (n = 0; n < res_n; ++n) {
        aln_res *p = res+n;
        for (i = 0; i < p->l_n; ++i) {
            if (p->la[i].merg_msg.x != 1) continue;
            // primary res
            for (j = 0; j <= p->la[i].cur_res_n; ++j) {
                all++;
                sam_flag = p->la[i].res[j].nsrand?(0x0):(0x10);
                // CIGAR
                if (prim_flag == 0 || AP.supp_soft) {
                    // QNAME/FLAG/RNAME/POS/MAPQ 
                    ksprintf(&sam_str, "%s\t%d\t%s\t%lld\t%d\t", name, sam_flag, bns->anns[p->la[i].res[j].chr-1].name, (long long)p->la[i].res[j].offset, p->la[i].mapQ); 
                    for (l = 0; l < p->la[i].res[j].cigar_len; ++l)
                        ksprintf(&sam_str, "%d%c", p->la[i].res[j].cigar[l]>>4, CIGAR_STR[p->la[i].res[j].cigar[l]&0xf]);
                } else {
                    sam_flag |= 0x800;
                    ksprintf(&sam_str, "%s\t%d\t%s\t%lld\t%d\t", name, sam_flag, bns->anns[p->la[i].res[j].chr-1].name, (long long)p->la[i].res[j].offset, p->la[i].mapQ); 
                    for (l = 0; l < p->la[i].res[j].cigar_len; ++l)
                        ksprintf(&sam_str, "%d%c", p->la[i].res[j].cigar[l]>>4, CIGAR_STR_HC[p->la[i].res[j].cigar[l]&0xf]);
                }
                // mate infomation
                kputs("\t*\t0\t0", &sam_str);
                if (prim_flag == 0 || AP.supp_soft) { // SOFT_CLIP
                    // SEQ and QUAL
                    int si;
                    if (p->la[i].res[j].nsrand == 1) {                  
                        kputc('\t', &sam_str);
                        kputs(seq, &sam_str);
                        kputc('\t', &sam_str);
                        if (qual && AP.comm) kputs(qual, &sam_str);
                        else kputc('*', &sam_str);
                    } else { // reverse strand
                        kputc('\t', &sam_str);
                        for (si = res->read_len-1; si >= 0; --si) kputc(com_nst_nt256_table[(uint8_t)seq[si]], &sam_str);
                        kputc('\t', &sam_str);
                        if (qual && AP.comm) for (si = res->read_len-1; si >=0; ++si) kputc(qual[si], &sam_str);   
                        else kputc('*', &sam_str);
                    }
                    prim_flag = 1;
                } else { // HARD_CLIP
                    // SEQ and QUAL
                    int si;
                    if (p->la[i].res[j].nsrand == 1) {                  
                        kputc('\t', &sam_str);
                        for (si = p->la[i].res[j].reg_beg-1; si < p->la[i].res[j].reg_end; ++si) kputc(seq[si], &sam_str);
                        kputc('\t', &sam_str);
                        if (qual && AP.comm) for (si = p->la[i].res[j].reg_beg-1; si < p->la[i].res[j].reg_end; ++si) kputc(qual[si], &sam_str);
                        else kputc('*', &sam_str);
                    } else { // reverse strand
                        kputc('\t', &sam_str);
                        for (si = p->la[i].res[j].reg_end-1; si >= p->la[i].res[j].reg_beg-1; --si) kputc(com_nst_nt256_table[(uint8_t)seq[si]], &sam_str);
                        kputc('\t', &sam_str);
                        if (qual && AP.comm) for (si = p->la[i].res[j].reg_end-1; si >= p->la[i].res[j].reg_beg-1; --si) kputc(qual[si], &sam_str);   
                        else kputc('*', &sam_str);
                    }
                }
                // optional_tags
                ksprintf(&sam_str, "\tNM:i:%d\tAS:i:%d", p->la[i].res[j].NM, p->la[i].res[j].score);
                check_cigar(p->la[i].res[j].cigar, p->la[i].res[j].cigar_len, name, res->read_len); 
                // alternative res XA
                if (j == 0) {
                    int m;
                    kstring_t XA_str; memset(&XA_str, 0, sizeof(kstring_t));
                    for (l = 0; l < p->la[i].XA_n; ++l) {
                        res_t *XA_res = p->la[i].XA[l];
                        ksprintf(&XA_str, "%s,%c%lld,", bns->anns[XA_res->chr-1].name,"-+"[XA_res->nsrand], (long long)XA_res->offset); 
                        for (m = 0; m < XA_res->cigar_len; m++)
                            ksprintf(&XA_str, "%d%c", (int)(XA_res->cigar[m]>>4), CIGAR_STR[(int)(XA_res->cigar[m] & 0xf)]);
                        check_cigar(XA_res->cigar, XA_res->cigar_len, name, res->read_len); 
                        ksprintf(&XA_str, ",%d;", XA_res->NM);
                    }
                    if (XA_str.l > 0) { 
                        ksprintf(&sam_str, "\tXA:Z:%s", XA_str.s);
                        free(XA_str.s); 
                    }
                }
                kputc('\n', &sam_str);
                fprintf(outp, "%s", sam_str.s);
                sam_str.l = 0;
            }
        }
    }
    if (all == 0) { // unmap
        ksprintf(&sam_str, "%s\t%d\t*\t*\t%d\t*\t*\t0\t0\t", name, 4, 0);
        kputs(seq, &sam_str);
        kputc('\t', &sam_str);
        if (qual) kputs(qual, &sam_str);
        else kputc('*', &sam_str);
        kputc('\n', &sam_str);
        fprintf(outp, "%s", sam_str.s);
    }
    free(sam_str.s);
}

void lamsa_aln_output(lamsa_aln_para AP, lamsa_seq_t *lamsa_seqs, kseq_t *seqs, int n_seqs, bntseq_t *bns)
{
    int i;
    for (i = 0; i < n_seqs; ++i) {
        lamsa_seq_t *la_seqs = lamsa_seqs+i;
        kseq_t *p = seqs+i;
        aln_res_output(AP, la_seqs->a_res, 3, p->name.s, p->seq.s, p->qual.s, bns);
    }
}

#ifdef __DEBUG__
#define CHUNK_SIZE 1
//#define CHUNK_READ_N 1
#endif
int lamsa_aln_core(const char *read_prefix, char *seed_result, seed_msg *s_msg, bwt_t *bwt, bntseq_t *bns, uint8_t *pac, lamsa_aln_para *AP)
{
    lamsa_seq_t *lamsa_seqs; int n_seqs, i;
    gzFile readfp; kseq_t *read_seq_t; FILE *seed_mapfp;
    // initialization for read.fa and seed.out, seed_msg
    if ((seed_mapfp = fopen(seed_result, "r")) == NULL) { fprintf(stderr, "\n[lamsa_aln_core] Can't open read file %s.\n", read_prefix); exit(1); }
    if ((readfp = gzopen(read_prefix, "r")) == NULL) { fprintf(stderr, "\n[lamsa_aln_core] Can't open seed-result file %s.\n", seed_result); exit(1); }
    kstream_t *fs = ks_init(readfp);
    read_seq_t = (kseq_t*)calloc(CHUNK_READ_N, sizeof(kseq_t));
    for (i = 0; i < CHUNK_READ_N; ++i) read_seq_t[i].f = fs;
    s_msg->read_count = 0;

    // alloc and initialization for auxiliary data
    thread_aux_t *aux;
    if (AP->n_thread < 1) AP->n_thread = 1;
    aux = (thread_aux_t*)calloc(AP->n_thread, sizeof(thread_aux_t));
    for (i = 0; i < AP->n_thread; ++i) {
        aux[i].tid = i; 
        aux[i].AP = AP; aux[i].bwt = bwt; aux[i].pac = pac; aux[i].bns = bns;
        aux_dp_init(aux+i, s_msg, *AP);
    }
    lamsa_seqs = lamsa_seq_init(CHUNK_READ_N, s_msg->seed_max, AP->res_mul_max);
    pthread_rwlock_init(&RWLOCK, NULL);

    // core loop
    char *gem_line = (char*)malloc(LINE_SIZE);
    while ((n_seqs = lamsa_read_seq(lamsa_seqs, read_seq_t, seed_mapfp, gem_line, LINE_SIZE, s_msg, CHUNK_READ_N)) != 0) { // read a chunk of read and other input data
        THREAD_READ_I = 0;
        if (AP->n_thread <= 1) {
            aux->n_seqs = n_seqs; 
            aux->w_seqs = read_seq_t;
            aux->lamsa_seqs = lamsa_seqs;
            lamsa_main_aln(aux);
        } else {
            pthread_t *tid; pthread_attr_t attr; 
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            tid = (pthread_t*)calloc(AP->n_thread, sizeof(pthread_t));
            int j;
            for (j = 0; j < AP->n_thread; ++j) {
                aux[j].n_seqs = n_seqs; 
                aux[j].w_seqs =read_seq_t;
                aux[j].lamsa_seqs = lamsa_seqs;
                pthread_create(&tid[j], &attr, lamsa_thread_aln, aux+j);
            }
            for (j = 0; j < AP->n_thread; ++j) pthread_join(tid[j], 0);
            free(tid);
        }
        // output align result 
        lamsa_aln_output(*AP, lamsa_seqs, read_seq_t, n_seqs, bns);
    }
    // free seqs
    lamsa_free_read_seq(lamsa_seqs, CHUNK_READ_N);
    pthread_rwlock_destroy(&RWLOCK);
    for (i = 0; i < AP->n_thread; ++i) aux_dp_free(aux+i, s_msg, AP);
    free(aux); gzclose(readfp); fclose(seed_mapfp); 
    free(read_seq_t); free(fs); free(gem_line);
    return 0;
}

int lamsa_gem(const char *ref_prefix, const char *read_prefix, char *bin_dir, int per_loci, int n_thread)
{
    char cmd[1024];
    sprintf(cmd, "bash %s/gem/gem_map.sh %s %s.seed %d %d", bin_dir, ref_prefix, read_prefix, per_loci, n_thread);
    fprintf(stderr, "[lamsa_aln] Executing gem-mapper ... \n");
    fprintf(stderr, "[lamsa_aln] Time consumption of gem-mapper:");
    if (system(cmd) != 0) { fprintf(stderr, "[lamsa_aln] Seeding undone, gem-mapper exit abnormally.\n"); exit(1); }
    fprintf(stderr, "[lamsa_aln] gem-mapper done!\n");
    return 0;
}

int lamsa_bwa(const char *ref_prefix, const char *read_prefix, char *bin_dir)
{
    char cmd[1024];
    sprintf(cmd, "bash %s/bwa/bwa_aln.sh %s %s.seed", bin_dir, ref_prefix, read_prefix);
    fprintf(stderr, "[lamsa_aln] Executing bwa aln ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lamsa_aln] Seeding undone, bwa aln exit abnormally.\n"); exit(1); }
    fprintf(stderr, "done!\n");
    return 0;
}

int lamsa_soap2dp(const char *ref_prefix, const char *read_prefix, char *bin_dir)
{
    char cmd[1024];
    sprintf(cmd, "bash %s/soap2-dp/soap2dp_aln.sh %s %s", bin_dir, ref_prefix, read_prefix);
    fprintf(stderr, "[lamsa_aln] Executing soap2dp aln ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lamsa_aln] Seeding undone, soap2dp aln exit abnormalloy.\n"); exit(1); }
    fprintf(stderr, "done!\n");
    return 0;
}

void print_sam_header(FILE *outp, const bntseq_t *bns)
{
	int i;
    for (i = 0; i < bns->n_seqs; ++i)
        fprintf(outp, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	fprintf(outp, "%s\n", lamsa_pg);
}

int lamsa_aln_c(const char *ref_prefix, const char *read_prefix, int seed_info, int seed_program, int no_seed_aln, char *seed_result, lamsa_aln_para *AP, char *bin_dir)
{
    clock_t t = clock();

    /* load index */
    fprintf(stderr, "[lamsa_aln] Restoring ref-indices ... ");
    bwt_t *bwt; bntseq_t *bns;
    char *str = (char*)calloc(strlen(ref_prefix)+5, 1);
    strcpy(str, ref_prefix); strcat(str, ".bwt"); bwt = bwt_restore_bwt(str);
    strcpy(str, ref_prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
    free(str);

    bns = bns_restore(ref_prefix);
    uint8_t *pac = (uint8_t*)calloc(bns->l_pac/4+1, 1);
    fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);	fprintf(stderr, "done!\n");


    /* split-seeding */
    seed_msg *s_msg;
    s_msg = seed_init_msg();
    if (seed_info) split_seed_info(read_prefix, *AP, s_msg);
    else split_seed(read_prefix, *AP, s_msg);

    if (!strcmp(seed_result, ""))
    {
        strcpy(seed_result, read_prefix);
        if (seed_program == 0)
            strcat(seed_result, ".seed.gem.map");
        else if (seed_program == 1)
            strcat(seed_result, ".seed.bwa.sam");
        else if (seed_program == 2)
            strcat(seed_result, ".seed.out.0");
        else { fprintf(stderr, "[lamsa_aln] Unknown seeding program option.\n"); return lamsa_aln_usage(); }
    }
    //excute soap2-dp program
    if (!no_seed_aln) 
    {
        if (seed_program == 0) lamsa_gem(ref_prefix, read_prefix, bin_dir, AP->per_aln_m, AP->n_thread);
        else if (seed_program == 1) lamsa_bwa(ref_prefix, read_prefix, bin_dir);
        else if (seed_program == 2) lamsa_soap2dp(ref_prefix, read_prefix, bin_dir);
        else { fprintf(stderr, "[lamsa_aln] Unknown seeding program option.\n"); return lamsa_aln_usage(); }
    }

    fprintf(stderr, "[lamsa_aln] Mapping reads to genome ...\n");
    print_sam_header(AP->outp, bns);

    lamsa_aln_core(read_prefix, seed_result, s_msg, bwt, bns, pac, AP);

    fprintf(stderr, "[lamsa_aln] Mapping done!\n");
    fprintf(stderr, "[lamsa_aln] Time consumption %.3f sec.\n", (float)(clock() -t )/CLOCKS_PER_SEC);

    seed_free_msg(s_msg);
    bwt_destroy(bwt); bns_destroy(bns); free(pac);
    return 0;
}

void init_aln_para(lamsa_aln_para *AP)
{
    AP->n_thread = 1;

    AP->seed_len = SEED_LEN;
    AP->seed_step = SEED_INTERVAL;
    //AP->seed_inv = SEED_INTERVAL;
    //AP->match_dis = AP->seed_step/25; // 4% indel allowed
    AP->per_aln_m = SEED_PER_LOCI;
    AP->first_loci_thd = SEED_FIRST_ROUND_LOCI;

    AP->SV_len_thd = SV_MAX_LEN;
    AP->ske_max = MAX_SKEL;
    AP->ovlp_rat = (float)OVLP_RAT;

    AP->split_len = SPLIT_ALN_LEN;
    AP->split_pen = SPLIT_ALN_PEN;
    AP->res_mul_max = RES_MAX_N;

    AP->hash_len = HASH_LEN;
    AP->hash_key_len = HASH_KEY;
    AP->hash_step = HASH_STEP;
    AP->hash_size = (int)pow(NT_N, AP->hash_key_len);
    
    AP->bwt_seed_len = BWT_KMER;
    AP->bwt_max_len = MAX_BWT_REG;

    AP->supp_soft = 0;
    AP->comm = 0;
    AP->outp = stdout;

    AP->frag_score_table = f_BCC_score_table;
    AP->match = MAT_SCORE; AP->mis = MIS_PEN;
    AP->gapo = OPEN_PEN; AP->gape = EXT_PEN;
    AP->end_bonus = 5;
    AP->zdrop = 100;
}

void lamsa_fill_mat(int mat, int mis, int8_t sc_mat[25])
{
    int i, j, k;
    for (i = k = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j)
            sc_mat[k++] = i==j ? mat : -mis;
        sc_mat[k++] = -1; // 'N'
    }
    for (j = 0; j < 5; ++j) sc_mat[k++] = -1; // 'N'
}

extern char *get_bin_dir(char *bin);

    //while ((c =getopt(argc, argv, "t:l:i:p:V:v:s:R:k:m:M:O:E:r:g:SCo:hH")) >= 0)
const struct option long_opt [] = {
    { "thread", 1, NULL, 't' },
    { "seed-len", 1, NULL, 'l' },
    { "seed-inv", 1, NULL, 'i' },
    { "max-loci", 1, NULL, 'p' },
    { "SV-len", 1, NULL, 'V' },
    { "ovlp-rat", 1, NULL, 'v' },
    { "max-skel", 1, NULL, 's' },
    { "max-reg", 1, NULL, 'R' },
    { "bwt-kmer", 1, NULL, 'k' },

    { "match-sc", 1, NULL, 'm' },
    { "mis-pen", 1, NULL, 'M' },
    { "open-pen", 1, NULL, 'O' },
    { "ext-pen", 1, NULL, 'E' },

    { "max-out", 1, NULL, 'r' },
    { "gap-split", 1, NULL, 'g' },
    { "soft-clip", 0, NULL, 'S' },
    { "comment", 0, NULL, 'C' },
    { "output", 1, NULL, 'o' },

    { "help", 0, NULL, 'h' },
    { "HELP", 0, NULL, 'H' },
    { 0, 0, 0, 0}
};

int lamsa_aln(int argc, char *argv[])
{
    int c;
    char *ref, *read;
    // parameters
    lamsa_aln_para *AP = (lamsa_aln_para*)calloc(1, sizeof(lamsa_aln_para));
    init_aln_para(AP);
    int no_seed_aln=0, seed_info=0, seed_program=0;
    char seed_result_f[1024]="";

    //while ((c =getopt(argc, argv, "t:m:M:O:E:S:r:V:g:s:l:i:p:d:o:NIA:")) >= 0)
    while ((c =getopt_long(argc, argv, "t:l:i:p:V:v:s:R:k:m:M:O:E:r:g:SCo:hHNI", long_opt, NULL)) >= 0)
    {
        switch (c)
        {
			case 't': AP->n_thread = atoi(optarg); break;

			case 'l': AP->seed_len = atoi(optarg); break;
			case 'i': AP->seed_step = atoi(optarg); break;
			case 'p': AP->per_aln_m = atoi(optarg); break;
            case 'V': AP->SV_len_thd = atoi(optarg); break;
            case 'v': if ((AP->ovlp_rat = atof(optarg)) < 0 || AP->ovlp_rat > 1) return lamsa_aln_usage(); break;
            case 's': AP->ske_max = atoi(optarg); break;
            case 'R': AP->bwt_max_len = atoi(optarg); break;
            case 'k': AP->bwt_seed_len = atoi(optarg); break;

            case 'm': AP->match = atoi(optarg); break;
            case 'M': AP->mis = atoi(optarg); break;
            case 'O': AP->gapo = atoi(optarg); break;
            case 'E': AP->gape = atoi(optarg); break;

            case 'r': AP->res_mul_max = atoi(optarg); break;
            case 'g': AP->split_len = atoi(optarg);
            case 'S': AP->supp_soft = 1; break;
            case 'C': AP->comm = 1; break;
            case 'o': AP->outp = fopen(optarg, "w"); 
                      if (AP->outp == NULL) { fprintf(stderr, "[lamsa_aln] Can not open output file: %s.\n", optarg); exit(1); }
                      break;
            case 'h': return lamsa_aln_usage();
            case 'H': return lamsa_aln_de_usage();
            case 'N': no_seed_aln = 1; break;
            case 'I': seed_info = 1; break; 
            default: return lamsa_aln_usage();
        }
    }
    AP->seed_inv = AP->seed_step-AP->seed_len;
    AP->match_dis = AP->seed_step / 25; // 4% indel allowed
    lamsa_fill_mat(AP->match, AP->mis, AP->sc_mat);
    if (argc - optind != 3)
        return lamsa_aln_usage();

    ref = strdup(argv[optind+1]);
    read =strdup(argv[optind+2]);

    char *bin_dir = get_bin_dir(argv[0]);

    lamsa_aln_c(ref, read, seed_info, seed_program, no_seed_aln, seed_result_f, AP, bin_dir);

    fclose(AP->outp); free(AP);
    free(ref); free(read);
    return 0;
}
