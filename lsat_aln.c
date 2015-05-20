#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include "lsat_aln.h"
#include "bwt.h"
#include "bntseq.h"
#include "frag_check.h"
#include "split_mapping.h"
#include "bwt_aln.h"
#include "lsat_heap.h"
#include "kseq.h"
#include "kstring.h"
#include "gem_parse.h"
#include "./lsat_sam_parse/sam.h"

KSEQ_INIT(gzFile, gzread)

int lsat_aln_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   lsat aln [options] <ref_prefix> <in.fa/fq>\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "         -t [INT]      The align type, <Best Coverage and Connect(1)> or <All Valid Result(2)>. [Def=1]\n\n");

    fprintf(stderr, "         -m [INT]      Match score for SW-alignment. [Def=1]\n");
    fprintf(stderr, "         -M [INT]      Mismatch penalty for SW-alignment. [Def=3]\n");
    fprintf(stderr, "         -O [INT]      Gap open penalty for SW-alignment. [Def=5]\n");
    fprintf(stderr, "         -E [INT]      Gap extension penalty for SW-alignment. [Def=2]\n");
    fprintf(stderr, "         -S [INT]      The score penalty of split-alignment. [Def=%d]\n\n", SPLIT_ALN_PEN);

    fprintf(stderr, "         -r [INT]      The maximun number of output results for a specific read region. [Def=%d]\n", RES_MAX_N);
    fprintf(stderr, "         -V [INT]      The expected maximun length of SV. [Def=%d]\n", SV_MAX_LEN);
    fprintf(stderr, "         -g [INT]      The minimum length of gap that causes a split-alignment. [Def=%d]\n\n", SPLIT_ALN_LEN);

    fprintf(stderr, "         -s [INT]      The seeding program, <gem(0)>, <bwa(1)> or <soap2-dp(2)>. [Def=0]\n");
	fprintf(stderr, "         -l [INT]      The seed length. [Def=Auto]\n");
	fprintf(stderr, "         -v [INT]      The interval length of adjacent seeds. [Def=Auto]\n");
	fprintf(stderr, "         -p [INT]      The maximun allowed number of a seed's locations. [Def=Auto]\n");
	fprintf(stderr, "         -d [INT]      The maximun number of seed's locations for first round's DP. [Def=Auto]\n\n");
	fprintf(stderr, "                         LSAT will set seed-len and seed-inv for every read based on their read-length.\n");
	fprintf(stderr, "                         Note: All of these 4 parameters(l,v,p,d) above should always be \"Auto\"\n");
	fprintf(stderr, "                               or set with specific value simultaneously.\n\n");
    
    fprintf(stderr, "         -o [STR]      The output file (SAM format). [Def=\"prefix_out.sam\"]\n\n");

    fprintf(stderr, "         -N            Do NOT excute seeding program, when seeds' alignment result existed already.\n");
    fprintf(stderr, "         -I            Seed information file has already existed.\n");
    fprintf(stderr, "         -A [STR]      The seeds' alignment result. When '-N' is used. [Def=\"seed_prefix.out.0\"]\n");
    fprintf(stderr, "\n");
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

int lsat_read_leve_n = 10;
int lsat_read_len[9]  = {2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};

int lsat_seed_len[10] = {  50,   50,   50,   75,   75,   75,  100,  100, 100, 100};

////                         40 31-46 30-40 40-50 40-48 40-47 47-54 46-52  45-50 40  
//int lsat_seed_inv[10] = {   0,   15,   50,   25,   50,   75,   50,   75, 100, 150}; 
//                         40 27-40 30-40 40-50 40-48 40-47 47-54 46-52  45-50 40  
int lsat_seed_inv[10] = {   0,   25,   50,   25,   50,   75,   50,   75, 100, 150}; 
int lsat_per_aln[10]  = { 200,  200,  200,  150,  150,  150,  100,  100, 100, 100}; //XXX
int lsat_mis_len[10]  = {  30,   30,   30,   25,   25,   25,   20,   20,  20,  20};//XXX
int lsat_min_thd[10]  = {   2,    2,    1,    1,    1,    1,    1,    1,   1,   1};

// for debug
char READ_NAME[1024];

int get_read_level(int read_len)
{
    int l;
    for (l = 0; l < lsat_read_leve_n-1; ++l)
        if (read_len <= lsat_read_len[l]) return l;
    return l;
}

seed_msg *seed_init_msg(void)
{
    seed_msg *msg = (seed_msg*)malloc(sizeof(seed_msg));

    msg->read_all = 0;
    msg->read_m = READ_MAX_NUM;

    msg->seed_len = (int *)calloc(READ_MAX_NUM, sizeof(int));
    msg->seed_inv = (int *)calloc(READ_MAX_NUM, sizeof(int));
    msg->seed_all = (int *)calloc(READ_MAX_NUM, sizeof(int));
    msg->read_name = (char **)malloc(READ_MAX_NUM * sizeof(char*));
    int i;
    for (i = 0; i < READ_MAX_NUM; ++i)
        msg->read_name[i] = (char*)malloc(1024 * sizeof(char));
    msg->last_len = (int *)malloc(READ_MAX_NUM * sizeof(int));
    msg->read_len = (int *)malloc(READ_MAX_NUM * sizeof(int));
    msg->read_level = (int*)malloc(READ_MAX_NUM * sizeof(int));
    msg->seed_max = 0;
    msg->read_max_len = 0;

    return msg;
}

void seed_free_msg(seed_msg *msg)
{
    int i;
    for (i = 0; i < msg->read_m; ++i) free(msg->read_name[i]);
    free(msg->read_name);
    free(msg->seed_len);
    free(msg->seed_inv);
    free(msg->seed_all);
    free(msg->last_len);
    free(msg->read_len);
    free(msg->read_level);
    free(msg);
}

int set_seed_argv(int *seed_len, int *seed_inv, int *read_level, int read_len)
{
    *read_level = get_read_level(read_len);
    *seed_len = lsat_seed_len[*read_level];
    *seed_inv = lsat_seed_inv[*read_level];
    return 0;
}

int split_seed(const char *prefix, seed_msg *s_msg)
{
    gzFile infp;
    kseq_t *seq;
    char out_f[1024], seed_head[1024], seed_seq[1024], seed_info[1024];
    FILE *outfp, *infofp;
    int m_read, seed_all, i;
    void *new_p;

    if ((infp = gzopen(prefix, "r")) == NULL) {
        fprintf(stderr, "[lsat_aln] Can't open read file %s\n", prefix);
        exit(-1);
    }
    seq = kseq_init(infp);

    strcpy(out_f, prefix); strcat(out_f, ".seed");
    if ((outfp = fopen(out_f, "w")) == NULL) {
        fprintf(stderr, "[lsat_aln] Can't open seed file %s\n", out_f);
        exit(-1);
    }
    strcpy(seed_info, prefix); strcat(seed_info, ".seed.info");
    if ((infofp = fopen(seed_info, "w")) == NULL) {
        fprintf(stderr, "[lsat_aln] Can't open seed info file %s\n", seed_info);
        exit(-1);
    }

    fprintf(stderr, "[lsat_aln] Spliting seed ... ");

    int seed_len, seed_inv, read_level;
    m_read = s_msg->read_m;
    while (kseq_read(seq) >= 0)
    {
        /* calculate seed-len and seed-interval based on read-length */
        set_seed_argv(&seed_len, &seed_inv, &read_level, seq->seq.l);
        seed_all = (1+ (seq->seq.l - seed_len) / (seed_len + seed_inv));
        seed_seq[seed_len] = '\n';
        if (seed_all > s_msg->seed_max) s_msg->seed_max = seed_all;
        if (seq->seq.l > s_msg->read_max_len) s_msg->read_max_len = seq->seq.l; //XXX
        if (s_msg->read_all == m_read-1)
        {
            m_read <<= 1;
            if ((new_p = (int*)realloc(s_msg->seed_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_len);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for seed_len[].\n");
                exit(1);
            }
            s_msg->seed_len = new_p;
            if ((new_p = (int*)realloc(s_msg->seed_inv , m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_inv);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for seed_inv[].\n");
                exit(-1);
            }
            s_msg->seed_inv = new_p;

            if ((new_p = (int*)realloc(s_msg->seed_all, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_all);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for seed_all[].\n");
                exit(1);
            }
            s_msg->seed_all = new_p;
            if ((new_p = (int*)realloc(s_msg->last_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->last_len);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for last_len[].\n");
                exit(-1);
            }
            s_msg->last_len = new_p;
            if ((new_p = (int*)realloc(s_msg->read_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_len);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_len = new_p;
            if ((new_p = (int*)realloc(s_msg->read_level, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_level);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for read_level[].\n");
                exit(-1);
            }
            s_msg->read_level = new_p;

            if ((new_p = (char**)realloc(s_msg->read_name, m_read * sizeof(char*))) == NULL)
            {
                free(s_msg->read_name);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_name = new_p;
            int i;
            for (i = m_read>>1; i < m_read; ++i)
                s_msg->read_name[i] = (char*)malloc(1024 * sizeof(char));
        }
        ++s_msg->read_all;
        s_msg->seed_len[s_msg->read_all] = seed_len;
        s_msg->seed_inv[s_msg->read_all] = seed_inv;
        s_msg->seed_all[s_msg->read_all] = s_msg->seed_all[s_msg->read_all-1] + seed_all;
        s_msg->last_len[s_msg->read_all] = seq->seq.l - seed_all * seed_len - (seed_all-1) * seed_inv;
        s_msg->read_len[s_msg->read_all] = seq->seq.l;
        s_msg->read_level[s_msg->read_all] = read_level;
        strcpy(s_msg->read_name[s_msg->read_all], seq->name.s);
        s_msg->read_m = m_read;

        for (i = 0; i < seed_all; ++i)
        {
            sprintf(seed_head, ">%s_%d:%d\n", seq->name.s, i, i*(seed_len+seed_inv));
            strncpy(seed_seq, seq->seq.s+i*(seed_len+seed_inv), seed_len);
            seed_seq[seed_len+1] = '\0';
            fputs(seed_head, outfp);
            fputs(seed_seq, outfp);
        }
        fprintf(infofp, "%s %d %d %d %d %d %d\n", seq->name.s, seed_len, seed_inv, read_level, seed_all, s_msg->last_len[s_msg->read_all], (int)seq->seq.l);
    }

    fprintf(stderr, "done.\n");
    gzclose(infp);
    fclose(outfp);
    fclose(infofp);
    kseq_destroy(seq);

    return 0;
}

int split_seed_info(const char *prefix, seed_msg *s_msg)
{
    char seed_info[1024];
    char read_name[1024];
    FILE *infofp;
    int m_read, seed_all, last_len, len, n, seed_len, seed_inv, read_level;
    void *new_p;

    strcpy(seed_info, prefix); strcat(seed_info, ".seed.info");
    if ((infofp = fopen(seed_info, "r")) == NULL)
    {
        fprintf(stderr, "[split seed] Can't open %s.\n", seed_info); 
        exit(-1);
    }
    m_read = s_msg->read_m;
    fprintf(stderr, "[lsat_aln] Parsing seeds' information ... ");
    
    while ((n = fscanf(infofp, "%s %d %d %d %d %d %d", read_name, &seed_len, &seed_inv, &read_level, &seed_all, &last_len, &len)) != EOF)
    {
        if (n != 7)
        {
            fprintf(stderr, "\n[split seed] INFO file error.[2]\n");
            exit(-1);
        }
        if (seed_all > s_msg->seed_max) s_msg->seed_max = seed_all;
        if (len > s_msg->read_max_len) s_msg->read_max_len = len;
        if (s_msg->read_all == m_read-1)
        {
            m_read <<= 1;
            if ((new_p = (int*)realloc(s_msg->seed_len , m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_len);
                fprintf(stderr, "\n[lsat aln] Can't allocate more memory for seed_len[].\n");
                exit(-1);
            }
            s_msg->seed_len =new_p;
            if ((new_p = (int*)realloc(s_msg->seed_inv , m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_inv);
                fprintf(stderr, "\n[lsat aln] Can't allocate more memory for seed_inv[].\n");
                exit(-1);
            }
            s_msg->seed_inv =new_p;
            if ((new_p = (int*)realloc(s_msg->seed_all, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->seed_all);
                fprintf(stderr, "\n[lsat aln] Can't allocate more memory for seed_all[].\n");
                exit(-1);
            }
            s_msg->seed_all =new_p;
            if ((new_p = (int*)realloc(s_msg->last_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->last_len);
                fprintf(stderr, "\n[lsat aln] Can't allocate more memory for last_len[].\n");
                exit(-1);
            }
            s_msg->last_len = new_p;
            if ((new_p = (int*)realloc(s_msg->read_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_len);
                fprintf(stderr, "\n[lsat aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_len = new_p;
            if ((new_p = (int*)realloc(s_msg->read_level, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_level);
                fprintf(stderr, "\n[lsat aln] Can't allocate more memory for read_level[].\n");
                exit(-1);
            }
            s_msg->read_level = new_p;
            //read_name
            if ((new_p = (char**)realloc(s_msg->read_name, m_read * sizeof(char*))) == NULL)
            {
                free(s_msg->read_name);
                fprintf(stderr, "\n[lsat_aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_name = new_p;
            int i;
            for (i = m_read>>1; i < m_read; ++i)
                s_msg->read_name[i] = (char*)malloc(1024 * sizeof(char));
        }
        ++s_msg->read_all;
        s_msg->seed_len[s_msg->read_all] = seed_len;
        s_msg->seed_inv[s_msg->read_all] = seed_inv;
        s_msg->seed_all[s_msg->read_all] = s_msg->seed_all[s_msg->read_all-1] + seed_all;
        s_msg->last_len[s_msg->read_all] = last_len;
        s_msg->read_len[s_msg->read_all] = len;
        s_msg->read_level[s_msg->read_all] = read_level;
        strcpy(s_msg->read_name[s_msg->read_all], read_name);
        s_msg->read_m = m_read;
        if (last_len != len - seed_all * seed_len - (seed_all-1)*seed_inv)
        {
            fprintf(stderr, "\n[split seed] INFO file error.[3]\n");
            exit(-1);
        }
    }

    fprintf(stderr, "done.\n");
    fclose(infofp);
    return 0;
}

aln_res *aln_init_res(int l_m, int read_len)
{
    int i, j;
    aln_res *a_res;
    a_res = (aln_res*)malloc(sizeof(aln_res));
    a_res->l_m = l_m, a_res->l_n = 0;
    a_res->la = (line_aln_res*)malloc(l_m * sizeof(line_aln_res));
    a_res->read_len = read_len;
    for (i = 0; i < l_m; ++i) {
        a_res->la[i].res_m = 10, a_res->la[i].cur_res_n = 0, a_res->la[i].split_flag = 0;
        a_res->la[i].res = (res_t*)malloc(10 * sizeof(res_t));
        for (j = 0; j < 10; ++j) {
            a_res->la[i].res[j].c_m = 100;
            a_res->la[i].res[j].cigar = (cigar32_t*)malloc(100 * sizeof(cigar32_t));
            a_res->la[i].res[j].cigar_len = 0;
        }
        a_res->la[i].tol_score = a_res->la[i].tol_NM = 0;
        a_res->la[i].trg_m = 10, a_res->la[i].trg_n = 0;
        a_res->la[i].trg = (line_node*)malloc(10 * sizeof(line_node));
        //a_res->la[i].merg_msg = ((*f_msg)+i)->merg_msg;
    }
    return a_res;
}

void aln_res_free(aln_res *res)
{
    int i,j,k;
    for (i = 0; i < res->l_m; ++i) {
        for (j = 0; j < res->la[i].res_m; ++j) {
            free(res->la[i].res[j].cigar);
        }
        free(res->la[i].res);
        free(res->la[i].trg);
    }
    free(res->la); free(res);
}

aln_reg *aln_init_reg(int read_len)
{
    aln_reg *reg = (aln_reg*)malloc(sizeof(aln_reg));
    reg->reg_n = 0; reg->reg_m = 1;
    reg->read_len = read_len;
    reg->reg = (reg_t*)malloc(sizeof(reg_t));
    return reg;
}

void aln_free_reg(aln_reg *reg) { free(reg->reg); free(reg); }

int reg_comp(const void *a, const void *b) { return (*(reg_t*)a).beg - (*(reg_t*)b).beg; }
void aln_sort_reg(aln_reg *a_reg) { qsort(a_reg->reg, a_reg->reg_n, sizeof(reg_t), reg_comp); }

void aln_merg_reg(aln_reg *a_reg, lsat_aln_para *AP) {
    //if (a_reg->reg_n == 0)return;
    int cur_i = 0, i;
    for (i = 1; i < a_reg->reg_n; ++i) {
		// merge
        if(a_reg->reg[i].beg - a_reg->reg[cur_i].end - 1 < 19) { //XXX
		    if(a_reg->reg[i].end > a_reg->reg[cur_i].end) 
				a_reg->reg[cur_i].end = a_reg->reg[i].end;
		}
        else {
            cur_i++;
            a_reg->reg[cur_i].beg = a_reg->reg[i].beg;
            a_reg->reg[cur_i].end = a_reg->reg[i].end;
        }
    }
    a_reg->reg_n = cur_i+1;
}

void push_reg(aln_reg *reg, reg_t r)
{
    if (reg->reg_n == reg->reg_m) {
        reg->reg_m <<= 1;
        if ((reg->reg = (reg_t*)realloc(reg->reg, reg->reg_m * sizeof(reg_t))) == NULL) {
            fprintf(stderr, "[push_reg] Not enough memory.\n");
            exit(0);
        }
    }
    reg->reg[reg->reg_n].beg = r.beg;
    reg->reg[reg->reg_n].end = r.end;
	reg->reg[reg->reg_n].refid = r.refid;
	reg->reg[reg->reg_n].is_rev = r.is_rev;
	reg->reg[reg->reg_n].ref_beg = r.ref_beg;
	reg->reg[reg->reg_n].ref_end = r.ref_end;
    reg->reg_n++;
}

// XXX dump remain-regions that are longger than L
int get_remain_reg(aln_reg *a_reg, aln_reg *remain_reg, lsat_aln_para *AP)
{
    if (a_reg->reg_n == 0) return 0;
    aln_sort_reg(a_reg); aln_merg_reg(a_reg, AP);
    int i, reg_thd=300; // XXX
    if (a_reg->reg[0].beg > AP->bwt_seed_len && a_reg->reg[0].beg-1 <= reg_thd)
		push_reg(remain_reg, (reg_t){a_reg->reg[0].refid, a_reg->reg[0].is_rev, 
				                     0, a_reg->reg[0].ref_beg, 
				                     1, a_reg->reg[0].beg-1});
    for (i = 1; i < a_reg->reg_n; ++i) {
        if (a_reg->reg[i].beg - a_reg->reg[i-1].end > AP->bwt_seed_len && a_reg->reg[i].beg-1-a_reg->reg[i-1].end <= reg_thd) {
			if (a_reg->reg[i].refid == a_reg->reg[i-1].refid && a_reg->reg[i].is_rev == a_reg->reg[i-1].is_rev)
				push_reg(remain_reg, (reg_t){a_reg->reg[i].refid, a_reg->reg[i].is_rev,
						                     a_reg->reg[i-1].ref_end, a_reg->reg[i].ref_beg, 
						                     a_reg->reg[i-1].end+1, a_reg->reg[i].beg-1});
			else 
				push_reg(remain_reg, (reg_t){-1, -1, 0, 0, a_reg->reg[i-1].end+1, a_reg->reg[i].beg-1});
		}
    }
	if (a_reg->read_len - a_reg->reg[i-1].end > AP->bwt_seed_len && a_reg->read_len-a_reg->reg[i-1].end <= reg_thd) 
		push_reg(remain_reg, (reg_t){a_reg->reg[i-1].refid, a_reg->reg[i-1].is_rev,
				                     a_reg->reg[i-1].ref_end, 0,
				                     a_reg->reg[i-1].end+1, a_reg->read_len});
    return remain_reg->reg_n;
}

void push_reg_res(aln_reg *reg, res_t res)
{
    if (reg->reg_n == reg->reg_m) {
        reg->reg_m <<= 1;
        if ((reg->reg = (reg_t*)realloc(reg->reg, reg->reg_m * sizeof(reg_t))) == NULL) {
            fprintf(stderr, "[push_reg_res] Not enough memory.\n");
            exit(0);
        }
    }
	reg->reg[reg->reg_n].refid = res.chr-1;
	reg->reg[reg->reg_n].is_rev = 1-res.nsrand;
    if (res.nsrand == 1) { // '+'
        reg->reg[reg->reg_n].beg = (res.cigar[0] & 0xf) == CSOFT_CLIP ? ((res.cigar[0]>>4)+1):1;
        reg->reg[reg->reg_n].end = (res.cigar[res.cigar_len-1] & 0xf) == CSOFT_CLIP ? (reg->read_len - (res.cigar[res.cigar_len-1]>>4)):reg->read_len;
		reg->reg[reg->reg_n].ref_beg = res.offset;
		reg->reg[reg->reg_n].ref_end = res.offset+refInCigar(res.cigar, res.cigar_len)-1;
    } else { // '-'
        reg->reg[reg->reg_n].beg = (res.cigar[res.cigar_len-1] & 0xf) == CSOFT_CLIP ? ((res.cigar[res.cigar_len-1]>>4) + 1):1;
        reg->reg[reg->reg_n].end = (res.cigar[0] & 0xf) == CSOFT_CLIP ? (reg->read_len - (res.cigar[0]>>4)):reg->read_len;
		reg->reg[reg->reg_n].ref_end = res.offset;
		reg->reg[reg->reg_n].ref_beg = res.offset+refInCigar(res.cigar, res.cigar_len)-1;
    }
    reg->reg_n++;
}

/*void aln_res_output(aln_res *res, aln_reg *reg, lsat_aln_per_para *APP) 
{
    int i, j;
    for (i = 0; i < res->l_n; ++i) {
        if (res->la[i].merg_msg.x==1) {
            for (j = 0; j <= res->la[i].cur_res_n; ++j) {
                // primary res
                fprintf(stdout, "%s\t%d\t%c\t%lld\t", APP->read_name, res->la[i].res[j].chr, "-++"[res->la[i].res[j].nsrand], (long long)res->la[i].res[j].offset); 
                printcigar(stdout, res->la[i].res[j].cigar, res->la[i].res[j].cigar_len); 
                fprintf(stdout, "\tAS:i:%d\tNM:i:%d", res->la[i].res[j].score, res->la[i].res[j].NM);
                if (reg) push_reg_res(reg, res->la[i].res[j]);
                check_cigar(res->la[i].res[j].cigar, res->la[i].res[j].cigar_len, APP->read_name, APP->read_len); 
                if (j == 0 && res->la[i].merg_msg.y != -1) {
                    int alt = res->la[i].merg_msg.y;
                    int k;
                    // alternative res
                    for (k = 0; k <= res->la[alt].cur_res_n; ++k) {
                        fprintf(stdout, "\tXA:Z:%d,%c%lld,", res->la[alt].res[k].chr, "-+"[res->la[alt].res[k].nsrand], (long long)res->la[alt].res[k].offset); 
                        printcigar(stdout, res->la[alt].res[k].cigar, res->la[alt].res[k].cigar_len); 
                        fprintf(stdout, ",%d;", res->la[alt].res[k].NM);
                    }
                }
                fprintf(stdout, "\n"); 
            }
        }
    }
}*/

//merg_msg: {1, 0} -> NOT merged or ONLY best
//          {1, 1} -> merged, best
//          {2, i} -> merged, alternative and best is `i`
//          {0,-1} -> dumped
void aln_res_output(aln_res *res, aln_reg *reg, lsat_aln_per_para *APP)
{
    int i, j;
    for (i = 0; i < res->l_n; ++i) {
        if (res->la[i].merg_msg.x == 1) {
            // primary res
			for (j = 0; j <= res->la[i].cur_res_n; ++j) {
				if (res->la[i].res[j].dump) continue;
				fprintf(stdout, "%s\t%d\t%c\t%lld\t", APP->read_name, res->la[i].res[j].chr, "-++"[res->la[i].res[j].nsrand], (long long)res->la[i].res[j].offset); 
				printcigar(stdout, res->la[i].res[j].cigar, res->la[i].res[j].cigar_len); 
				fprintf(stdout, "\tAS:i:%d\tNM:i:%d", res->la[i].res[j].score, res->la[i].res[j].NM);
				if (reg) push_reg_res(reg, res->la[i].res[j]);
				check_cigar(res->la[i].res[j].cigar, res->la[i].res[j].cigar_len, APP->read_name, APP->read_len); 
				// alternative res
				if (j == 0 && res->la[i].merg_msg.y == 1) {
					int k,l;
					for (l = 0; l < res->l_n; ++l) {
						if (res->la[l].merg_msg.x == 2 && res->la[l].merg_msg.y == i) {
							for (k = 0; k <= res->la[l].cur_res_n; ++k) {
								if (res->la[l].res[j].dump) continue;
								fprintf(stdout, "\tXA:Z:%d,%c%lld,", res->la[l].res[k].chr, "-+"[res->la[l].res[k].nsrand], (long long)res->la[l].res[k].offset); 
								printcigar(stdout, res->la[l].res[k].cigar, res->la[l].res[k].cigar_len); 
								fprintf(stdout, ",%d;", res->la[l].res[k].NM);
							}
						}
					}
				}
				fprintf(stdout, "\n"); 
			}
        }
    }
}

sam_msg *sam_init_msg(void)
{
    sam_msg *s = (sam_msg*)malloc(sizeof(sam_msg));
    s->sam_n = 0;
    s->sam_m = 1;
    s->sam = (sam_t*)malloc(sizeof(sam_t));
    s->sam->cigar_s = (kstring_t*)calloc(1, sizeof(kstring_t));
    s->sam->NM = 0;
    return s;
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
            msg[i].at[j].cigar = (cigar32_t*)malloc(7 * sizeof(cigar32_t));//XXX default value for 3-ed
            msg[i].at[j].cigar_len = 0;
            msg[i].at[j].cmax = 7;
            msg[i].at[j].bmax = 0;
        }
    }
    return msg;
}

void aln_realc_msg(aln_msg *a_msg, int seed_max, int per_aln_M, int per_aln_m)
{
    int i,j;
    for (i = 0; i < seed_max; ++i)
    {
        a_msg[i].at = (aln_t*)realloc(a_msg[i].at, per_aln_m * sizeof(aln_t));
        if (a_msg[i].at == NULL)
        {
            fprintf(stderr, "\n[lsat_aln] Not enough memory.\n");
            exit(0);
        }
        for (j = per_aln_M; j < per_aln_m; ++j)
        {
            a_msg[i].at[j].cigar = (cigar32_t*)malloc(7 * sizeof(cigar32_t));//XXX default value for 3-ed
            a_msg[i].at[j].cigar_len = 0;
            a_msg[i].at[j].cmax = 7;
            a_msg[i].at[j].bmax = 0;
        }
    }
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
            (*f_node)[i][j].son_max = 100; // son_max XXX
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

void fnode_realc(frag_dp_node ***f_node, int seed_m, int per_aln_M, int per_aln_m)
{
    int i, j;
    for (i = 0; i < seed_m; ++i)
    {
        (*f_node)[i] = (frag_dp_node*)realloc((*f_node)[i], per_aln_m * sizeof(frag_dp_node));
        for (j = per_aln_M; j < per_aln_m; ++j) {
            (*f_node)[i][j].seed_i = i, (*f_node)[i][j].aln_i = j;
            (*f_node)[i][j].son_max = 100; // son_max XXX
            (*f_node)[i][j].son = (line_node*)calloc(100, sizeof(line_node));
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
            fprintf(stderr, "\n[lsat_aln] Cigar ERROR 1.\n");
            fprintf(stderr, "%s\n",s);
            exit(-1);
        }
        op = toupper(*t);
        switch (op)
        {
            case 'M':	op = CMATCH;	break;
            case 'I':	op = CINS;		bi += x;	break;
            case 'D':	op = CDEL;		bd += x;	break;
                        //case 'S':	op = CSOFT_CLIP;		bi += x;	break;
                        //case 'H':	op = CHARD_CLIP;		bi += x;	break;
            default:	fprintf(stderr, "\n[lsat_aln] Cigar ERROR 2.\n"); exit(-1); break;
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

void setAmsg(aln_msg *a_msg, int32_t read_x, int aln_y, 
        int read_id, sam_t sam)
    //   char srand, int NM, char *cigar)
{   //read_x: (除去unmap和repeat的)read序号, aln_y: read对应的比对结果序号(从1开始)
    a_msg[read_x-1].read_id = read_id;			//from 1
    a_msg[read_x-1].at[aln_y-1].chr = sam.chr;
    a_msg[read_x-1].at[aln_y-1].offset = sam.offset;	//1-base
    a_msg[read_x-1].at[aln_y-1].nsrand = ((sam.nsrand=='+')?1:-1);
    a_msg[read_x-1].at[aln_y-1].NM = sam.NM;
    a_msg[read_x-1].n_aln = aln_y;
    setCigar(a_msg, read_x-1, aln_y-1,  sam.cigar_s->s);
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
void set_aln_msg(aln_msg *a_msg, int32_t read_x, int aln_y, int read_id, map_t map)
{   
    a_msg[read_x-1].read_id = read_id;			//from 1
    a_msg[read_x-1].at[aln_y-1].chr = (map.chr[3] == 'X' ? 23 : (map.chr[3] == 'Y' ? 24 :atoi(map.chr+3))); //XXX
    a_msg[read_x-1].at[aln_y-1].offset = map.offset;	//1-base
    a_msg[read_x-1].at[aln_y-1].nsrand = ((map.strand=='+')?1:-1);
    a_msg[read_x-1].at[aln_y-1].NM = map.NM;
    a_msg[read_x-1].n_aln = aln_y;
    set_cigar(a_msg[read_x-1].at+aln_y-1, map.cigar);
}

void init_aln_per_para(lsat_aln_per_para *APP, seed_msg *s_msg, int read_n)
{
    strcpy(APP->read_name, s_msg->read_name[read_n]);
    APP->read_len = s_msg->read_len[read_n];

    APP->seed_len = s_msg->seed_len[read_n];
    APP->seed_inv = s_msg->seed_inv[read_n];
    APP->read_level = s_msg->read_level[read_n];
    APP->seed_step = APP->seed_len + APP->seed_inv;
    APP->last_len = s_msg->last_len[read_n];

    APP->seed_all = s_msg->seed_all[read_n] - s_msg->seed_all[read_n-1];

    APP->per_aln_n = lsat_per_aln[APP->read_level]; // depend on seed_len
    APP->min_thd = lsat_min_thd[APP->read_level]; // depend on seed_len XXX
    APP->frag_score_table = f_BCC_score_table;

    APP->match_dis = APP->seed_step/25; // depend on seed_len, seed_inv, 4% indel allowed XXX
}

void init_aln_para(lsat_aln_para *AP)
{
    AP->aln_type = 1;   // BCC
    AP->per_aln_m = PER_ALN_N; 

    AP->SV_len_thd = SV_MAX_LEN;
    AP->split_len = SPLIT_ALN_LEN;
    AP->split_pen = SPLIT_ALN_PEN;

    AP->res_mul_max = RES_MAX_N;

    AP->hash_len = HASH_LEN;
    AP->hash_key_len = HASH_KEY;
    AP->hash_step = HASH_STEP;
    
    AP->bwt_seed_len = 19; //XXX

    AP->gapo = 5; AP->gape = 2;
    AP->match = 1; AP->mis = 3;
}

//for best coverage and connect
//add "overlap ins"
int get_fseed_dis(aln_msg *a_msg, int pre, int pre_a, int i, int j, int *flag, lsat_aln_per_para *APP, lsat_aln_para *AP)    //(i,j)对应节点，来自pre的第pre_a个aln
{
    if (pre == -1 || i == -1) { *flag = F_MATCH; return 0; }	//for bound node
    if (pre == i) {// XXX overlap???
        if (pre_a == j) *flag = F_MATCH; 
        else *flag = F_UNCONNECT;
        return 0;
    }
    if (a_msg[i].at[j].chr != a_msg[pre].at[pre_a].chr || a_msg[i].at[j].nsrand != a_msg[pre].at[pre_a].nsrand)	//different chr or different srnad
    {
        *flag = F_CHR_DIF; return 0;//PRICE_DIF_CHR;
    }

    int seed_len = APP->seed_len; int seed_inv = APP->seed_inv;
    int64_t exp = a_msg[pre].at[pre_a].offset + a_msg[pre].at[pre_a].nsrand * (a_msg[i].read_id - a_msg[pre].read_id) * (seed_len+seed_inv);	
    int64_t act = a_msg[i].at[j].offset;
    int64_t dis = a_msg[pre].at[pre_a].nsrand * ((a_msg[pre].read_id < a_msg[i].read_id)?(act-exp):(exp-act)) - (((a_msg[pre].at[pre_a].nsrand) * (a_msg[pre].read_id-a_msg[i].read_id) < 0)?(a_msg[pre].at[pre_a].len_dif):(a_msg[i].at[j].len_dif));

    int mat_dis = APP->match_dis;
    if (dis <= mat_dis && dis >= -mat_dis) {
        if (abs(a_msg[pre].read_id - a_msg[i].read_id) == 1) *flag = F_MATCH;
        else if (abs(a_msg[pre].read_id - a_msg[i].read_id) < 10) *flag = F_MISMATCH; // XXX long dis
        else *flag = F_LONG_MISMATCH;
    } else if (dis > mat_dis && dis < AP->SV_len_thd) *flag = F_DELETE;
    //else if (dis < -mat_dis && dis >= (0-(abs(a_msg[i].read_id-a_msg[pre].read_id)*(seed_len+seed_inv)-seed_len))) *flag = F_INSERT; // nonoverlap
    else if (dis < -mat_dis && dis >= -AP->SV_len_thd) *flag = F_INSERT; // overlap ins XXX
    else { *flag = F_UNCONNECT; return 0; }
    dis=abs(dis); dis += (((a_msg[pre].at[pre_a].nsrand) * (a_msg[pre].read_id - a_msg[i].read_id) < 0)? a_msg[i].at[j].cigar_len : a_msg[pre].at[pre_a].cigar_len);
    return dis; 
}

void fnode_set(frag_dp_node *f_node, line_node from, 
               int score, int dis_pen, 
               uint8_t match_flag, int dp_flag)
{
    (*f_node).son_flag = F_INIT;
    (*f_node).from = from;
    (*f_node).score = score;
    (*f_node).dis_pen = dis_pen;

    (*f_node).match_flag = match_flag;
    (*f_node).dp_flag = dp_flag;

    // general init
    (*f_node).node_n = 1;
    (*f_node).in_de = 0;
    (*f_node).son_n = 0;
    (*f_node).max_score = score;
    (*f_node).max_node.x = (*f_node).seed_i, (*f_node).max_node.y = (*f_node).aln_i;
    (*f_node).trg_n = 0;
}

//node_score = sum(node_num) + sum(connect_score) XXX
int frag_dp_init(frag_dp_node **f_node, 
                 aln_msg *a_msg, int seed_i, 
                 line_node from, 
                 lsat_aln_per_para *APP, lsat_aln_para *AP,
                 int dp_flag)
{
    int i;
    if (from.x == START_NODE.x) { //UNLIMITED
        for (i = 0; i < a_msg[seed_i].n_aln; ++i)
            //init-score: 0 == 1 + 1 - F_SV_PEN
            fnode_set(&(f_node[seed_i][i]), from, 1, a_msg[seed_i].at[i].cigar_len, F_MATCH, dp_flag);
    } else {
        int con_flag, con_score, dis;
        for (i = 0; i < a_msg[seed_i].n_aln; ++i) {
            
            dis = get_fseed_dis(a_msg, from.x, from.y, seed_i, i, &con_flag, APP, AP);
            //											XXX MATCH, MISMATCH: same score?
            //											mismatch too long, socre??? XXX
            con_score = APP->frag_score_table[con_flag];

            if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
                fnode_set(&(f_node[seed_i][i]), from, 2+con_score, dis, con_flag, dp_flag);
            else f_node[seed_i][i].dp_flag = 0-dp_flag;
        }
    }
    return 0;
}

int frag_dp_per_init(frag_dp_node **f_node, aln_msg *a_msg, 
                     int seed_i, int aln_i, line_node from, 
                     lsat_aln_per_para *APP, lsat_aln_para *AP,
                     int dp_flag)
{
    if (from.x == START_NODE.x) { //UNLIMITED
        //init-score: 0 == 1 + 1 - F_SV_PEN
        fnode_set(&(f_node[seed_i][aln_i]), from, 1, a_msg[seed_i].at[aln_i].cigar_len, F_MATCH, dp_flag);
    } else {
        int con_flag, con_score, dis;
        dis = get_fseed_dis(a_msg, from.x, from.y, seed_i, aln_i, &con_flag, APP, AP);
        //											XXX MATCH, MISMATCH: same score?
        //											mismatch too long, socre??? XXX
        con_score = APP->frag_score_table[con_flag];

        if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
            fnode_set(&(f_node[seed_i][aln_i]), from, 2+con_score, dis, con_flag, dp_flag);
        else f_node[seed_i][aln_i].dp_flag = 0-dp_flag;
    }
    return 0;
}

int fnode_add_son(frag_dp_node **f_node,
                  line_node fa, line_node son)
{
    ++f_node[fa.x][fa.y].in_de;
    int son_n = f_node[fa.x][fa.y].son_n;
    if (son_n == f_node[fa.x][fa.y].son_max) {
        f_node[fa.x][fa.y].son_max <<= 1;
        f_node[fa.x][fa.y].son = (line_node*)realloc(f_node[fa.x][fa.y].son, f_node[fa.x][fa.y].son_max * sizeof(line_node));
    }
    f_node[fa.x][fa.y].son[son_n] = son;
    ++f_node[fa.x][fa.y].son_n;
    return 0;
}

//pruning	XXX
//for best coverage and connect
//    similar to MEM, match/mismatch have the highest priority in dp.
int frag_dp_update(frag_dp_node **f_node, aln_msg *a_msg, 
                   int seed_i, int aln_i, int start, 
                   lsat_aln_per_para *APP, lsat_aln_para *AP, 
                   int dp_flag)
{
    int i, j, con_flag, con_score, dis;
    line_node max_from;
    int max_score, max_flag, max_dis;

    max_from = f_node[seed_i][aln_i].from;
    max_score = f_node[seed_i][aln_i].score;
    max_dis = 0; max_flag = 0;

    for (i = seed_i - 1; i >= start; --i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == dp_flag)// || f_node[i][j].dp_flag == UPDATE_FLAG)
            {
                dis = get_fseed_dis(a_msg, i, j, seed_i, aln_i, &con_flag, APP, AP);
                if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
                    continue;
                // '+' strand: successor-Match
                //if (con_flag > f_node[i][j].son_flag+1) //
                 //   continue;
                if (a_msg[i].at[j].nsrand == 1) {
                    if (f_node[i][j].son_flag <= F_MISMATCH) {
                        continue;
                        /// cigar-diff
                        /*if (con_flag > F_MISMATCH)
                            continue;
                        else { 
                            max_from = (line_node){i,j};
                            con_score = APP->frag_score_table[con_flag];
                            max_score = f_node[i][j].score + 1 + con_score;
                            max_flag = con_flag;
                            max_dis = f_node[i][j].dis_pen + dis;
                            goto UPDATE;
                        }*/
                    } 
                }
                // '-' strand: precursor-Match 
                ///XXX LONG_MISMATCH
                ///XXX cigar-diff, when score equal with each other
                else if (con_flag <= F_MISMATCH) {
                    max_from = (line_node){i,j};
                    con_score = APP->frag_score_table[con_flag];
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                    goto UPDATE;
                }

                con_score = APP->frag_score_table[con_flag];
                if (f_node[i][j].score + 1 + con_score > max_score)	{
                    max_from = (line_node){i,j};
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                } else if (f_node[i][j].score + 1 + con_score == max_score) {
                    if (f_node[i][j].dis_pen + dis < max_dis)
                    {
                        max_from = (line_node){i,j};
                        max_score = f_node[i][j].score + 1 + con_score;
                        max_flag = con_flag;
                        max_dis = f_node[i][j].dis_pen + dis;
                    }
                }
            }
        }
    }
UPDATE:
    if (max_from.x != f_node[seed_i][aln_i].from.x || max_from.y != f_node[seed_i][aln_i].from.y) {
        f_node[max_from.x][max_from.y].son_flag = max_flag;
        f_node[seed_i][aln_i].from = max_from;
        f_node[seed_i][aln_i].score = max_score;
        f_node[seed_i][aln_i].dis_pen = max_dis;
        f_node[seed_i][aln_i].match_flag = max_flag;
        f_node[seed_i][aln_i].node_n = f_node[max_from.x][max_from.y].node_n+1;

        fnode_add_son(f_node, max_from, (line_node){seed_i, aln_i}); //generate tree
    }
    return 0;
}

node_score *node_init_score(int n) {
    node_score *ns = (node_score*)malloc(sizeof(node_score));
    ns->max_n = n;
    ns->node_n = 0;
    ns->node = (line_node*)malloc(n * sizeof(line_node));
    ns->score = (int*)malloc(n * sizeof(int));
    return ns;
}

void node_free_score(node_score *ns) {
    free(ns->score); free(ns->node);
    free(ns);
}

int heap_add_node(node_score *ns, line_node node, int score)
{
    int i;
    if (ns->node_n < ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->node[ns->node_n++] = node;
        return -1;
    } else if (ns->node_n == ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->node[ns->node_n++] = node;
        build_node_min_heap(ns);
        return -1;
    } else { // (ns->node_n == max_n, min-heap is already built.
        return node_heap_update_min(ns, node, score); 
    }
}


//XXX 在输出node时就已经排序？即在dp过程中排序？
void line_sort_endpos(line_node **line, int *line_end,
                      trig_node **t_node, int *tri_n,
                      int li, int len) {
    int i, j, k, end, t;
    line_node tmp; trig_node tr;
    for (i = 0; i < len-1; ++i) {
        for (j = i+1; j < len; ++j) {
            if (line[j][line_end[j]-1].x > line[i][line_end[i]-1].x) {
                // swap_line
                end = line_end[j] > line_end[i] ? line_end[j] : line_end[i];
                for (k = 0; k < end+L_EXTRA; ++k) {
                    tmp = line[i][k], line[i][k] = line[j][k], line[j][k] = tmp;   
                }
                t = line_end[i], line_end[i] = line_end[j], line_end[j] = t;
                end = tri_n[j] > tri_n[i] ? tri_n[j] : tri_n[i];
                for (k = 0; k < end; ++k) {
                    tr = t_node[i][k], t_node[i][k] = t_node[j][k], t_node[j][k] = tr;   
                }
                t = tri_n[i], tri_n[i] = tri_n[j], tri_n[j] = t;
            }
        }
    }
}
void line_sort_endpos1(line_node **line, int *line_end,
                       int li, int len) {
    int i, j, k, end, t;
    line_node tmp;
    for (i = 0; i < len-1; ++i) {
        for (j = i+1; j < len; ++j) {
            if (line[j][line_end[j]-1].x > line[i][line_end[i]-1].x) {
                // swap_line
                end = line_end[j] > line_end[i] ? line_end[j] : line_end[i];
                for (k = 0; k < end+L_EXTRA; ++k) {
                    tmp = line[i][k], line[i][k] = line[j][k], line[j][k] = tmp;   
                }
                t = line_end[i], line_end[i] = line_end[j], line_end[j] = t;
            }
        }
    }
}

int line_merge(int a, int b, line_node **line, int *line_end) {
    int s1, e1, s2, e2, s, e, hi;
    float rat1, rat2;

    s2 = line[a][0].x, e2 = line[a][line_end[a]-1].x;
    ///if (line[b][line_end[b]+1].x == -1) { // NO merged-line
    if (L_MF(line,line_end,b) & L_NMERG) {// Not merged
        hi = b;
        s1 = line[b][0].x, e1 = line[b][line_end[b]-1].x;
    } else { // [b] was merged already
        hi = L_MH(line, line_end, b);
        s1 = L_LB(line, line_end, hi), e1 = L_RB(line, line_end, hi);
    }
    s = (s2 > s1) ? s2 : s1;
    e = (e2 < e1) ? e2 : e1;
    // new strategy: no mattter how big/small the overlap is, do merge.
    if (s <= e) {
        L_LB(line, line_end, hi) = s1+s2-s, L_RB(line, line_end, hi) = e1+e2-e;
        L_MF(line, line_end, hi) = L_MERGH;
        L_MF(line, line_end, a) = L_MERGB, L_MH(line, line_end, a) = hi;
        return 1;
    } else {
        L_MF(line, line_end, a) = L_NMERG;
        return 0;
    }
    /*
    rat1 = (e-s+1+0.0)/(e1-s1+0.0);
    rat2 = (e-s+1+0.0)/(e2-s2+0.0);
    if (rat1>=0.7 || rat2>=0.7) { //0.7? XXX
        ///line[hi][line_end[hi]+1].x = s1+s2-s, line[hi][line_end[hi]+1].y = e1+e2-e;
        L_LB(line, line_end, hi) = s1+s2-s, L_RB(line, line_end, hi) = e1+e2-e;
        L_MF(line, line_end, hi) = L_MERGH;
        ///line[a][line_end[a]+1].x = -2, line[a][line_end[a]+1].y = hi;
        L_MF(line, line_end, a) = L_MERGB, L_MH(line, line_end, a) = hi;
        return 1;
    }
    else {
        ///line[a][line_end[a]+1].x = -1;
        L_MF(line, line_end, a) = L_NMERG;
        return 0;
    }*/
}

// select best/secondary merged-line, based on line-node number OR based on line-score? XXX
// rule: 1. Better Score
//       2. Short line, if Score is equal
//            
// Best & Secondary: 1. Best is unique, 
//                   2. Secondary is allowed to be multi, 
//            		 3. Secondary is allowed to equal the Best
// new
// Best and all secondary whose score exceed best/2
void line_filter(line_node **line, int *line_end, int li, int len, trig_node **tri_node, int *tri_n, frag_dp_node **f_node, aln_msg *a_msg, int per_max_multi)
{
    int i, ii, j, k, l, m_ii;
    line_node **m_b = (line_node**)malloc(len * sizeof(line_node*)); // (i, score)
	int **m_f = (int**)malloc(len * sizeof(int*)); // merged line (after filtering)
	for (i = 0; i < len; ++i) {
		m_b[i] = (line_node*)malloc(len * sizeof(line_node));
		m_f[i] = (int*)malloc((len+1) * sizeof(int));
	}
	int *m_bn = (int*)malloc(len * sizeof(int)); 
	int *m_fn = (int*)malloc(len * sizeof(int));
    int m_i = -1, d_i=0;
    int b_score, s_score;
	int m_head, f_merge, min_l, max_r, head;
    // check merged-lines
    for (i = li; i < li+len; ++i) {
        if (L_MF(line, line_end, i) & L_NMERG) { // NOT merged
			++m_i;
			m_b[m_i][0].x = i, m_b[m_i][0].y = -2;
			m_bn[m_i] = 1;
		} else if (L_MF(line, line_end, i) & L_MERGH) { // merged, head
			++m_i;
            m_b[m_i][0].x = i, m_b[m_i][0].y = L_LS(line, line_end, i); // NO secondary
			m_bn[m_i] = 1;
		} else { // merged, body : inter-line candidate
            if (m_i < 0) { fprintf(stderr, "[line_filter] %s BUG.\n", READ_NAME); exit(-1); }
			m_b[m_i][m_bn[m_i]].x = i, m_b[m_i][m_bn[m_i]].y = L_LS(line, line_end, i);
			m_bn[m_i]++;
        }
    }
	// select best/secondary line, set flag, set inter-line
    node_score *ns = node_init_score(per_max_multi);
	for (i = 0; i <= m_i; ++i) {
		if (m_b[i][0].y == -2) { // not merged
			m_f[i][0] = m_b[i][0].x;
			m_fn[i] = 1;
			continue;
		}
		b_score = s_score = 0;
		for (j = 0; j < m_bn[i]; ++j) {
			if (m_b[i][j].y > b_score) {
				s_score = b_score;
				b_score = m_b[i][j].y;
			} else if (m_b[i][j].y > s_score)
				s_score = m_b[i][j].y;
		}
		if (s_score > b_score/2) f_merge = 1; // get M largest scores
		else f_merge = 0;                     // get the largest score

		m_fn[i] = 1; // m_f [0: best, 1: head, 1...: body
		if (f_merge) {
            //1. init_heap
            ns->node_n = 0;
			for (j = 0; j < m_bn[i]; ++j) {
                //5. output_M_nodes
				if (m_b[i][j].y > b_score/2) {
                    //2. add_node & build_heap & update_heap
                    int ret = heap_add_node(ns, (line_node){m_b[i][j].x,-1}, m_b[i][j].y);
                    if (ret == -2) { // dumped
                        L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
                        tri_n[m_b[i][j].x] = 0;
                    } else if (ret != -1) { // ret is dumped
                        L_MF(line, line_end, ret) |= L_DUMP;
                        tri_n[ret] = 0;
                    }
                } else {
					L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
					tri_n[m_b[i][j].x] = 0;
				}
			} 	
            //3. output_M_nodes, set head/body flags
            build_node_minpos_heap(ns);
            // head
            line_node head_n = node_heap_extract_minpos(ns);
            m_head = head_n.x;
            L_MF(line, line_end, m_head) = L_MERGH;

            min_l = line[m_head][0].x;
            max_r = line[m_head][line_end[m_head]-1].x;
            if (L_LS(line, line_end, m_head) == b_score) m_f[i][0] = m_head; // best
            m_f[i][m_fn[i]++] = m_head;

            int body;
            while ((body = node_heap_extract_minpos(ns).x) != -1) {
                L_MF(line, line_end, body) = L_MERGB;
                L_MH(line, line_end, body) = m_head;
                min_l = MINOFTWO(min_l, line[body][0].x);
                max_r = MAXOFTWO(max_r, line[body][line_end[body]-1].x);
                if (L_LS(line, line_end, body) == b_score) m_f[i][0] = body; // best
                m_f[i][m_fn[i]++] = body;
            }
            // set left/right boundary
			L_LB(line, line_end, m_head) = MINOFTWO(line[m_head][0].x, min_l);
			L_RB(line, line_end, m_head) = MAXOFTWO(line[m_head][line_end[m_head]-1].x, max_r);
		} else {
			for (j = 0; j < m_bn[i]; ++j) {
				if (m_b[i][j].y == b_score) {
					L_MF(line, line_end, m_b[i][j].x) =  L_NMERG;
					m_f[i][0] = m_b[i][j].x; // best
					m_f[i][m_fn[i]++] = m_b[i][j].x;
				} else {
					L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
					tri_n[m_b[i][j].x] = 0;
				}
			}
		}
		for (ii = 1; ii < m_fn[i]; ++ii) {
			j = m_f[i][ii];
			for (k = 0; k < tri_n[j]; ++k) {
                head = -1;
				// check inter-line: from j+1 -> merge_end
                for (l = j+1; l < li+len && (L_MF(line, line_end, l) & 0x3)==0; ++l) { // L_MERGB
                    if (line[l][0].x > tri_node[j][k].n1.x && line[l][line_end[l]-1].x < tri_node[j][k].n2.x) {
                        // check for MIS(INV)
						if (f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].match_flag == F_MISMATCH ||
                                f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].match_flag == F_LONG_MISMATCH) {
                            int x_s = line[l][0].x, y_s = line[l][0].y;
                            int x_e = line[l][line_end[l]-1].x, y_e = line[l][line_end[l]-1].y;
                            int x1 = tri_node[j][k].n1.x, y1 = tri_node[j][k].n1.y, x2 = tri_node[j][k].n2.x, y2 = tri_node[j][k].n2.y;
                            int nsrand = a_msg[x_s].at[y_s].nsrand;
                            if (nsrand == a_msg[x1].at[y1].nsrand ||
                                    a_msg[x_s].at[y_s].chr != a_msg[x1].at[y1].chr ||
                                    nsrand * a_msg[x_s].at[y_s].offset < nsrand * a_msg[x2].at[y2].offset ||
                                    nsrand * a_msg[x_e].at[y_e].offset > nsrand * a_msg[x1].at[y1].offset)
                                continue;
                        }
                        E_LB(line, line_end, l) = tri_node[j][k].n1.x;
                        E_RB(line, line_end, l) = tri_node[j][k].n2.x;
                        L_MF(line, line_end, l) = L_INTER;
                        if (tri_node[j][k].tri_flag == 0x2 && (f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].trg_n & 0x2)) // pre
                            f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].trg_n -= 0x2;
                        else if (tri_node[j][k].tri_flag == 0x1 && (f_node[tri_node[j][k].n1.x][tri_node[j][k].n1.y].trg_n & 0x1)) // next
                            f_node[tri_node[j][k].n1.x][tri_node[j][k].n1.y].trg_n -= 0x1; 
                        if (head == -1) {
                            L_MF(line, line_end, l) |= L_NMERG;
                            head = l;
                        } else {
                            L_MF(line, line_end, l) |= L_MERGB;
                            L_MH(line, line_end, l) = head;
                            L_MF(line, line_end, head) = L_INTER|L_MERGH;
                        }
                    }
                }
            }
		}
	}
    node_free_score(ns);

    // remove candidate "trigger-line", left-most or right-most line
    if (m_i > 0) {
        int L;
        // left-boundary tri-line
        m_ii = 0;
		// best line in merged-lines
        l = line[m_f[m_ii][0]][line_end[m_f[m_ii][0]]-1].x - line[m_f[m_ii][0]][0].x;
        L = line[m_f[m_ii+1][0]][line_end[m_f[m_ii+1][0]]-1].x - line[m_f[m_ii+1][0]][0].x;
        if (l < 2 && L >= 2) { // DUMP
			for (i = 1; i < m_fn[m_ii]; ++i)
				L_MF(line, line_end, m_f[m_ii][i]) = L_DUMP;
        }
        // right-boundary tri-line
        m_ii = m_i;
        l = line[m_f[m_ii][0]][line_end[m_f[m_ii][0]]-1].x - line[m_f[m_ii][0]][0].x;
        L = line[m_f[m_ii-1][0]][line_end[m_f[m_ii-1][0]]-1].x - line[m_f[m_ii-1][0]][0].x;
        if (l < 2 && L >= 2) { // DUMP
			for (i = 1; i < m_fn[m_ii]; ++i)
				L_MF(line, line_end, m_f[m_ii][i]) = L_DUMP;
        }
    }
	for (i = 0; i < len; ++i) { free(m_b[i]); free(m_f[i]); }
    free(m_b); free(m_bn); free(m_f); free(m_fn);
}

void line_filter1(line_node **line, int *line_end, int li, int len, int per_max_multi)
{
    int i,j, m_ii;
	int *m_n = (int*)malloc(len * sizeof(int));
    line_node **m_b = (line_node**)malloc(len * sizeof(line_node*)); // (i, score)
	for (i = 0; i < len; ++i) m_b[i] = (line_node*)malloc(len * sizeof(line_node));
    int m_i = -1, d_i=0;
    int b_score, s_score;
    for (i = li; i < li+len; ++i) {
        if (L_MF(line, line_end, i) & L_NMERG) continue; // NOT merged
        else if (L_MF(line, line_end, i) & L_MERGH) {    // merged, head
            ++m_i;
            m_b[m_i][0].x = i, m_b[m_i][0].y = L_LS(line, line_end, i);
			m_n[m_i] = 1;
        } else { // merged, body
            if (m_i < 0) { fprintf(stderr, "[line_filter] %s BUG.\n", READ_NAME); exit(-1); }
			m_b[m_i][m_n[m_i]].x = i, m_b[m_i][m_n[m_i]].y = L_LS(line, line_end, i);
			m_n[m_i]++;
        }
    }
    int min_l, max_r, f_merge, m_head;
    node_score *ns = node_init_score(per_max_multi);
	for (i = 0; i <= m_i; ++i) { // m_b[i][0 ... m_n[i]-1]
		// get best & secondary score
		b_score = s_score = 0;
		for (j = 0; j < m_n[i]; ++j) {
			if (m_b[i][j].y > b_score) {
				s_score = b_score;
				b_score = m_b[i][j].y;
			} else if (m_b[i][j].y > s_score)
				s_score = m_b[i][j].y;
		}
		if (s_score > b_score/2) f_merge = 1;
		else f_merge = 0;
		// select best & secondary line, set flag
		if (f_merge) {
            // 1. init heap
            ns->node_n = 0;
			for (j = 0; j < m_n[i]; ++j) {
				if (m_b[i][j].y > b_score/2) {
                    // 2. add node & build heap & update heap
                    int ret = heap_add_node(ns, (line_node){m_b[i][j].x, -1}, m_b[i][j].y);
                    if (ret == -2) { // dumped
                        L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
                    } else if (ret != -1) { // ret is dumped
                        L_MF(line, line_end, ret) |= L_DUMP;
                    }
                } else L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
			} 	
            // 3. ouput_M_nodes, set head/body flags
            build_node_minpos_heap(ns);
            line_node head_n = node_heap_extract_minpos(ns);
            m_head = head_n.x;
            L_MF(line, line_end, m_head) = L_MERGH;

            min_l = line[m_head][0].x;
            max_r = line[m_head][line_end[m_head]-1].x;

            int body;
            while ((body = node_heap_extract_minpos(ns).x) != -1) {
                L_MF(line, line_end, body) = L_MERGB;
                L_MH(line, line_end, body) = m_head;
                min_l = MINOFTWO(min_l, line[body][0].x);
                max_r = MAXOFTWO(max_r, line[body][line_end[body]-1].x);
            }
            // set left/right boundary
			L_LB(line, line_end, m_head) = MINOFTWO(line[m_head][0].x, min_l);
			L_RB(line, line_end, m_head) = MAXOFTWO(line[m_head][line_end[m_head]-1].x, max_r);
		} else {
			for (j = 0; j < m_n[i]; ++j) {
				if (m_b[i][j].y == b_score)
					L_MF(line, line_end, m_b[i][j].x) =  L_NMERG;
				else L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
			}
		}
	}
    node_free_score(ns);
	for (i = 0; i < len; ++i) free(m_b[i]);
	free(m_b); free(m_n);
}

// internal mini merge ? XXX
//     same to the process in "set-bound"?
// internal-merge, based on inter-trigger
//   INS/INV
//   XXX check match_flag
//      MIS(INV): chr-dif and offset
void line_inter_bound(line_node **line, int *line_end, int li, int len, trig_node **tri_node, int *tri_n, frag_dp_node **f_node, aln_msg *a_msg)
{
    int i,j,head,k;
    for (j = li; j < li+len; ++j) {
        for (k = 0; k < tri_n[j]; ++k) {
            head = -1;
            for (i = li; i < li+len; ++i) {
                ///if (line[i][line_end[i]+1].x != -3) continue;
                if (!(L_MF(line, line_end, i) & L_DUMP)) continue;
                if (line[i][0].x > tri_node[j][k].n1.x && line[i][line_end[i]-1].x < tri_node[j][k].n2.x) {
                    // check for MIS(INV)
                    if (f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].match_flag == F_MISMATCH ||
                        f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].match_flag == F_LONG_MISMATCH) {
                        int x_s = line[i][0].x, y_s = line[i][0].y;
                        int x_e = line[i][line_end[i]-1].x, y_e = line[i][line_end[i]-1].y;
                        int x1 = tri_node[j][k].n1.x, y1 = tri_node[j][k].n1.y, x2 = tri_node[j][k].n2.x, y2 = tri_node[j][k].n2.y;
                        int nsrand = a_msg[x_s].at[y_s].nsrand;
                        if (nsrand == a_msg[x1].at[y1].nsrand ||
                            a_msg[x_s].at[y_s].chr != a_msg[x1].at[y1].chr ||
                            nsrand * a_msg[x_s].at[y_s].offset < nsrand * a_msg[x2].at[y2].offset ||
                            nsrand * a_msg[x_e].at[y_e].offset > nsrand * a_msg[x1].at[y1].offset)
                            continue;
                    }
                    ///line[i][line_end[i]].x = tri_node[j][k].x1;
                    ///line[i][line_end[i]].y = tri_node[j][k].x2;
                    E_LB(line, line_end, i) = tri_node[j][k].n1.x;
                    E_RB(line, line_end, i) = tri_node[j][k].n2.x;
                    L_MF(line, line_end, i) = L_INTER;
                    if (tri_node[j][k].tri_flag == 0x2 && (f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].trg_n & 0x2)) // pre
                        f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].trg_n -= 0x2;
                    else if (tri_node[j][k].tri_flag == 0x1 && (f_node[tri_node[j][k].n1.x][tri_node[j][k].n1.y].trg_n & 0x1)) // next
                        f_node[tri_node[j][k].n1.x][tri_node[j][k].n1.y].trg_n -= 0x1; 
                    if (head == -1) {
                        ///line[i][line_end[i]+1].x = -1;
                        L_MF(line, line_end, i) |= L_NMERG;
                        head = i;
                    } else {
                        ///line[i][line_end[i]+1] = (line_node){-2, head};
                        L_MF(line, line_end, i) |= L_MERGB;
                        L_MH(line, line_end, i) = head;
                        ///line[head][line_end[head]+1].x = 0;
                        L_MF(line, line_end, head) = L_INTER|L_MERGH;
                    }
                }
            }
        }
    }
}

int line_remove(line_node **line, int *line_end, int li, int len)
{
    int *head = (int*)malloc(len * sizeof(int));
    int l, j, cur_i=li;
    // remove line with "DUMP" flag
    for (l = li; l < li+len; ++l) {
        if (!(L_MF(line, line_end, l) & L_DUMP)) {
            if (L_MF(line, line_end, l) & L_MERGH)
                head[l-li] = cur_i;
            if (l == cur_i) {
                ++cur_i;
                continue;
            }
            line_end[cur_i] = line_end[l];
            for (j = 0; j < line_end[l]+L_EXTRA; ++j)
                line[cur_i][j] = line[l][j];
            //merged, body
            if (!(L_MF(line, line_end, cur_i)&L_NMERG) && !(L_MF(line, line_end, cur_i)&L_MERGH))
                    L_MH(line, line_end, cur_i) = head[L_MH(line, line_end, cur_i)-li];
            ++cur_i;
        }
    }
    free(head);
    return cur_i-li;
}

//XXX 将NotMerged与Merged统一起来
/* line: [0 ~ end-1] -> nodea; 
         [end] -> (left-b, right-b); 
         [end+1] -> (merged-flag, merge-head)/(start, end) 
         // (-1,):NOT merged
         // (-2,i):Merged and head is i
         // (-2,t):inter-merged, head is t
         // (-3,x):Merged, but dumped
         // (a,b): head, start from a to b
         [end+2] -> (line-score, x)
         */
int line_set_bound(line_node **line, int *line_end, 
                   int li, int *o_len, int left, int right, 
                   trig_node **t_node, int *tri_n, frag_dp_node **f_node, aln_msg *a_msg, int per_max_multi) 
{
    if (*o_len <= 0) return 0;
    int i, j, len = *o_len;
    // sort by end-pos
    line_sort_endpos(line, line_end, t_node, tri_n, li, len);
    // merge lines
    ///line[li][line_end[li]+1].x = -1;
    L_MF(line, line_end, li) = L_NMERG;
    for (i = 1; i < len; ++i)
        line_merge(li+i, li+i-1, line, line_end);
    // filter best/secondary line, and their tri-nodes, then find inter-line 
    line_filter(line, line_end, li, len, t_node, tri_n, f_node, a_msg, per_max_multi);
    *o_len = line_remove(line, line_end, li, len);
    len = *o_len;
    // set line-boundary
    int l,r,m,s,e;
    r = right;
    for (i = li; i < li+len; ++i) {
        if (!(L_MF(line, line_end, i) & L_DUMP)) {
            E_RB(line, line_end, i) = r;
            break;
        }
    }
    if (i == li+len) {fprintf(stderr, "[line_set_bound] %s BUG 0.\n", READ_NAME); exit(-1);}
    // set right-b
    if (L_MF(line, line_end, i) & L_NMERG)
        s = line[i][0].x, e = line[i][line_end[i]-1].x;
    else
        s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
    m = 0; //Merged head OR NOT-merged head
    for (++i; i < li+len; ++i) {
        if (L_MF(line, line_end, i) & (L_DUMP|L_INTER)) continue;
        if (L_MF(line, line_end, i) & L_NMERG) {
            r = s, l = line[i][line_end[i]-1].x;
            //set left-b for last-line(s)
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = line[i][0].x, e = line[i][line_end[i]-1].x;
            m = i;
        } else if (L_MF(line, line_end, i) & L_MERGH) {
            r = s, l = L_RB(line, line_end, i);
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
            m = i;
        } //merged and body
        E_RB(line, line_end, i) = r;
    }
    for (j = i-1; j >= m; --j) {
        if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
        E_LB(line, line_end, j) = left;
    }
    // set inter-bound
    //line_inter_bound(line, line_end, li, len, t_node, tri_n, f_node, a_msg);
    return 0;
}

int line_set_bound1(line_node **line, int *line_end, 
                   int li, int *o_len, int left, int right, int per_max_multi) 
{
    if (*o_len <= 0) return 0;
    int i, j, len=*o_len;
    // sort by end-pos
    line_sort_endpos1(line, line_end, li, len);
    // merge lines
    ///line[li][line_end[li]+1].x = -1;
    L_MF(line, line_end, li) = L_NMERG;
    for (i = 1; i < len; ++i)
        line_merge(li+i, li+i-1, line, line_end);
    // filter best/secondary line, and their tri-nodes, then find inter-line 
    line_filter1(line, line_end, li, len, per_max_multi);
    *o_len = line_remove(line, line_end, li, len);
    len = *o_len;
    // set line-boundary
    int l,r,m,s,e;
    r = right;
    for (i = li; i < li+len; ++i) {
        if (!(L_MF(line, line_end, i) & L_DUMP)) {
            E_RB(line, line_end, i) = r;
            break;
        }
    }
    if (i == li+len) {fprintf(stderr, "[line_set_bound] %s BUG 0.\n", READ_NAME); exit(-1);}
    // set right-b
    if (L_MF(line, line_end, i) & L_NMERG)
        s = line[i][0].x, e = line[i][line_end[i]-1].x;
    else
        s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
    m = 0; //Merged head OR NOT-merged head
    for (++i; i < li+len; ++i) {
        if (L_MF(line, line_end, i) & (L_DUMP|L_INTER)) continue;
        if (L_MF(line, line_end, i) & L_NMERG) {
            r = s, l = line[i][line_end[i]-1].x;
            //set left-b for last-line(s)
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = line[i][0].x, e = line[i][line_end[i]-1].x;
            m = i;
        } else if (L_MF(line, line_end, i) & L_MERGH) {
            r = s, l = L_RB(line, line_end, i);
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
            m = i;
        } //merged and body
        E_RB(line, line_end, i) = r;
    }
    for (j = i-1; j >= m; --j) {
        if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
        E_LB(line, line_end, j) = left;
    }
    // set inter-bound
    //line_inter_bound(line, line_end, li, len, t_node, tri_n, f_node, a_msg);
    return 0;
}

//TODO: save all (node,score)? OR just save the first `M` big (node,score)
void node_add_score(int score, line_node node, node_score *ns)
{
    // XXX
    //if (score <= 2) return;
    int i;
    if (ns->node_n <= ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->node[(ns->node_n)++] = node;
    } else {
        fprintf(stderr, "[lsat_aln] node_add_score ERROR. (%d %d)\n", ns->node_n, ns->max_n);
        exit(0);
    }
}

//1.F_MATCH/F_MISMATCH/F_LONG_MISMATCH
//2.F_INS/F_DEL:
//  a.Better-Score
//  b.Short-Distance
line_node get_max_son(frag_dp_node **f_node, int x, int y)
{
    line_node max, son;
    int i, max_score = 0, max_flag, max_dis=f_node[x][y].son[0].x - x;
    for (i = 0; i < f_node[x][y].son_n; ++i) {
        son = f_node[x][y].son[i];
        max_flag = f_node[son.x][son.y].match_flag;
        /// XXX LONG_MIS_MATCH
        if (max_flag == F_MATCH || max_flag == F_MISMATCH)// || max_flag == F_LONG_MISMATCH)
            return son;
        if (f_node[son.x][son.y].max_score > max_score || (f_node[son.x][son.y].max_score == max_score && son.x - x < max_dis)) {
            max = son;
            max_score = f_node[son.x][son.y].max_score;
            max_dis = son.x - x;
        } 
        /*if (son.x - x == max_dis) { // same dis, better score
            if (f_node[son.x][son.y].max_score > max_score)
            {
                max = son;
                max_score = f_node[son.x][son.y].max_score;
            }
        } else if (son.x - x < max_dis) {
            max = son;
            max_score = f_node[son.x][son.y].max_score;
            max_dis = son.x - x;
        }*/
    }
    return max;
}

void cut_branch(frag_dp_node **f_node, int x, int y, node_score *ns) {
    line_node max_son = get_max_son(f_node, x, y);
    line_node son, max_node;
    int i;
    for (i = 0; i < f_node[x][y].son_n; ++i) {
        son = f_node[x][y].son[i];
        if (nodeEq(max_son, son)) continue;
        //update score: 
        //New_Max = Old_Max - branch_score - con_score 
        //        = Old_Max - son_score + 1
        //Because: son_score = branch_score + 1 + con_score
        f_node[son.x][son.y].from = START_NODE;
        f_node[son.x][son.y].max_score -= (f_node[son.x][son.y].score - 1);
        max_node = f_node[son.x][son.y].max_node;
        f_node[max_node.x][max_node.y].node_n -= (f_node[son.x][son.y].node_n-1);
        node_add_score(f_node[son.x][son.y].max_score, max_node , ns);
    }
    // refresh max_socre/node of (x,y), based on max_son, consider negative edge
    if (f_node[x][y].score > f_node[max_son.x][max_son.y].max_score) { // negative edge
        f_node[max_son.x][max_son.y].in_de = -1;
        f_node[max_son.x][max_son.y].from = START_NODE;
        f_node[max_son.x][max_son.y].max_score -= (f_node[max_son.x][max_son.y].score - 1);
        max_node = f_node[max_son.x][max_son.y].max_node;
        f_node[max_node.x][max_node.y].node_n -= (f_node[max_son.x][max_son.y].node_n-1);
        node_add_score(f_node[max_son.x][max_son.y].max_score, max_node, ns);
        f_node[x][y].son_n = 0;
        f_node[x][y].max_node = (line_node){x,y};
        f_node[x][y].max_score = f_node[x][y].score;
    } else {
        f_node[x][y].son_n = 1;
        f_node[x][y].son[0] = max_son;
        f_node[x][y].max_node = f_node[max_son.x][max_son.y].max_node;
        f_node[x][y].max_score = f_node[max_son.x][max_son.y].max_score;
    }
    f_node[x][y].in_de = 0; // leaf-node, avaible for branch-track
}

//trackback from (x,y)
void branch_track_new(frag_dp_node **f_node, int x, int y, node_score *ns)
{
    int max_score;line_node max_node; 
    f_node[x][y].in_de = -1;
    if (f_node[x][y].son_n == 0) {
        max_node = f_node[x][y].max_node = (line_node){x,y};
        max_score = f_node[x][y].max_score = f_node[x][y].score;
    } else {
        max_node = f_node[x][y].max_node;
        max_score = f_node[x][y].max_score;
    }

    line_node fa = f_node[x][y].from, son; int fa_son_n;

    while (fa.x != START_NODE.x) {
        fa_son_n = f_node[fa.x][fa.y].son_n;
        if (fa_son_n == 1) { // 继续回溯
            if (f_node[fa.x][fa.y].score > max_score) { //negative edge
                son = f_node[fa.x][fa.y].son[0];
                f_node[son.x][son.y].in_de = -1;
                f_node[son.x][son.y].from = START_NODE;
                f_node[son.x][son.y].max_score -= (f_node[son.x][son.y].score - 1);
                f_node[max_node.x][max_node.y].node_n -= (f_node[son.x][son.y].node_n-1);
                node_add_score(f_node[son.x][son.y].max_score, max_node, ns);
                f_node[fa.x][fa.y].son_n = 0;
                max_score = f_node[fa.x][fa.y].score;
                max_node = fa;
            }
            f_node[fa.x][fa.y].max_score = max_score;
            f_node[fa.x][fa.y].max_node = max_node;
            f_node[fa.x][fa.y].in_de = -1;
            fa = f_node[fa.x][fa.y].from;
        } else { // fa_son_n > 1, TREE_BRANCH
            --f_node[fa.x][fa.y].in_de;
            if (f_node[fa.x][fa.y].in_de == 0)
                cut_branch(f_node, fa.x, fa.y, ns);
            return;
        }
    }
    // TREE_START; 
    node_add_score(max_score, max_node, ns);// add (max_node, max_score);
    return;
}
//trackback from (i,j)
/*line_node branch_track(frag_dp_node **f_node, int x, int y, node_score *ns)
{
    f_node[x][y].in_de = -1;
    line_node max_node = f_node[x][y].max_node = (line_node){x,y};
    int max_score = f_node[x][y].max_score = f_node[x][y].score;

    if (f_node[x][y].from.x == START_NODE.x)
        return max_node;

    line_node fa_node = f_node[x][y].from;
    int fa_son_n = f_node[fa_node.x][fa_node.y].son_n;

    while (fa_son_n == 1)
    {
        if (f_node[fa_node.x][fa_node.y].score > max_score) { //neg edge
            node_add_score(max_score, max_node, ns);
            max_score = f_node[fa_node.x][fa_node.y].score;
            max_node = fa_node;
        }
        f_node[fa_node.x][fa_node.y].max_score = max_score;
        f_node[fa_node.x][fa_node.y].max_node = max_node;
        f_node[fa_node.x][fa_node.y].in_de = -1;
        if (f_node[fa_node.x][fa_node.y].from.x == START_NODE.x)
            return fa_node;
        fa_node = f_node[fa_node.x][fa_node.y].from;
        fa_son_n = f_node[fa_node.x][fa_node.y].son_n;
    }
    --f_node[fa_node.x][fa_node.y].in_de;
    return fa_node;
}*/

//for multi-dp-lines
//line, line_end: start at 1
//left and right are 'MIN_FLAG'	XXX
int frag_mini_dp_multi_line(frag_dp_node **f_node, aln_msg *a_msg, 
                            lsat_aln_per_para *APP, lsat_aln_para *AP,
                            int left_b, int right_b, 
                            line_node **line, int *line_end)
{
    if (left_b+2 >= right_b) return 0;

    line_node head = START_NODE; // unlimited
    int start, end;
    int i, j, k, l, mini_dp_flag = WHOLE_FLAG;// (whole DP)
    int l_i, max_score, node_i;
    line_node left, right;
    line_node _right;
    //int candi_n; line_node candi[max_multi];

    left.x = left_b, right.x = right_b;
    start = left.x+1;
    end = right.x-1;

    //first dp int
    for (i = start; i <= end; ++i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            //if (f_node[i][j].dp_flag == mini_dp_flag || f_node[i][j].dp_flag == 0-mini_dp_flag)
                frag_dp_per_init(f_node, a_msg, i, j, head, APP, AP, mini_dp_flag);
        }
    }
    //dp update
    for (i = start+1; i <= end; ++i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag)
                frag_dp_update(f_node, a_msg, i, j, start, APP, AP, mini_dp_flag);
        }
    }
    // tree-pruning
    node_score *ns = node_init_score(AP->per_aln_m); //size of node_score XXX
    //node_score *ns = node_init_score(80); //size of node_score
    line_node fa_node;
    for (i = end; i >=start; --i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag && f_node[i][j].in_de == 0) //leaf node
                branch_track_new(f_node, i, j, ns); 
        }
    }
    l_i = 0;
    int score;
    while (1) {
        _right = node_pop(ns, &score);
        if (_right.x == START_NODE.x) break;
        node_i = f_node[_right.x][_right.y].node_n-1;
        line_end[l_i] = node_i+1;
        L_LS(line,line_end,l_i) = score;
        while (_right.x != head.x) {
            if (node_i < 0) { 
                fprintf(stderr, "\n[frag mini dp multi] %s node_i BUG 1.\n", APP->read_name); 
                exit(0); }
            line[l_i][node_i--] = _right;
            f_node[_right.x][_right.y].dp_flag = MULTI_OUT;
            _right = f_node[_right.x][_right.y].from;
        }
        if (node_i >= 0) { 
            fprintf(stderr, "\n[frag mini dp multi] %s node_i BUG 2.\n", APP->read_name); 
            exit(0); }
        ++l_i;
    }
    node_free_score(ns);
    return l_i;
    /*
    //
    // XXX tree-pruning?
    //store and sort candidate node, the biggest 'max-multi' cand-nodes
    candi_n=0;
    for (i = right.x-1; i > left.x; --i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            //candi-nodes
            if (f_node[i][j].backtrack_flag == DP_BACK_PEAK && f_node[i][j].dp_flag == mini_dp_flag && f_node[i][j].score > 0)
            {
                //sort candi-nodes in 'candi', big->small
                for (k = 0; k < candi_n; ++k) {
                    if (f_node[i][j].score > f_node[candi[k].x][candi[k].y].score)
                        break;
                }
                if (k >= max_multi) continue;
                if (candi_n < max_multi) candi_n++;
                for (l = candi_n-1; l > k; --l) {
                    candi[l] = candi[l-1];
                }
                candi[k].x = i; candi[k].y = j;
            }
        }
    }
    if (candi_n == 0) return 0;
    //find 'anchor-match' line
    //store in line[0], line_end[0]
    line_end[0] = 0;
    //backtrack
    //Uni-Best  or  All-Best
    l_i = 1;
    for (i = 0; i < candi_n; ++i) {
        if (i > 0) { //check for the backtrack_flag
            _right = candi[i];
            while (_right.x != head.x) {
                if (f_node[_right.x][0].backtrack_flag == DP_BACK_LOCATED)
                    goto NextLoop;
                _left = f_node[_right.x][_right.y].from;
                _right = _left;
            }
        }
        _right = candi[i];
        node_i = f_node[_right.x][_right.y].node_n - 1;
        line_end[l_i] = node_i+1;
        while (_right.x != head.x) {
            f_node[_right.x][0].backtrack_flag = DP_BACK_LOCATED;		//XXX
            if (node_i < 0) { fprintf(stderr, "\n[frag mini multi dp] node_i < 0, BUG.\n"); exit(0); }
            line[l_i][node_i] = _right;
            --node_i;
            _left = f_node[_right.x][_right.y].from;
            _right = _left;
        }
        if (node_i >= 0) { fprintf(stderr, "\n[frag mini multi dp] node_i >= 0, BUG.\n"); exit(0); }
        l_i++;
        NextLoop:;
    }*/
    /*for (i = 1; i < l_i; ++i) {
        fprintf(stdout, "mini line# %d:\n", i);
        for (j = 0; j < line_end[i]; ++j)
            fprintf(stdout, "(%d,%d)\t", line[i][j].x, line[i][j].y);
        printf("\n");
    }*/
    return l_i;
}

int trg_dp_line(frag_dp_node **f_node, aln_msg *a_msg, frag_msg **f_msg,
                lsat_aln_per_para *APP, lsat_aln_para *AP, 
                int left, int right, 
                line_node **line, int *line_end, 
                line_node **_line, int *_line_end, int per_max_multi) 
{
    int i, j;
    int l = frag_mini_dp_multi_line(f_node, a_msg, APP, AP, left, right, line, line_end);//, 0, 0);
    line_set_bound1(line, line_end, 0,  &l, left, right, per_max_multi);
    ///l = line_remove(line, line_end, 0, l);
    return l;
}
//mini-dp 过程中，种子比对结果被多重输出的问题
//解决办法:每次mini-dp结束后，修改该次输出的node对应的dp-flag，
//         以保证在后面的mini-dp中，它们不再被使用到.
//for best coverage and connect
//fill-up for inter-gaps in min-line: XXX 如何填补? 如何避免出现错误填补,即将不该填补的补上了
int frag_mini_dp_line(frag_dp_node **f_node, aln_msg *a_msg, 
                      lsat_aln_per_para *APP, lsat_aln_para *AP,
                      line_node left, line_node right, 
                      line_node *line, int *de_score,
                      int _head, int _tail)
{
    line_node head = ((_head)?left:START_NODE);
    int old_score;
    if (_head==0 || _tail==0) old_score = 1;
    else old_score = 2+f_BCC_score_table[f_node[right.x][right.y].match_flag];
    int i, j, mini_dp_flag = MULTI_FLAG;
    //dp init
    for (i = left.x + 1; i < right.x; ++i) {
        //mini_dp XXX?
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag || f_node[i][j].dp_flag == 0-mini_dp_flag)
                frag_dp_per_init(f_node, a_msg, i, j, head, APP, AP, mini_dp_flag);
        }
    }
    //dp update
    for (i = left.x + 2; i < right.x; ++i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag)
                frag_dp_update(f_node, a_msg, i, j, left.x+1, APP, AP, mini_dp_flag);
        }
    }
    //find backtrack start node
    int max_score, max_dis = 0, max_n = 0;
    line_node max_node = head;//left;
    //UNLIMITED 
    if (_tail == 0) {
        max_score = old_score;
        for (i = right.x - 1; i > left.x; --i) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                if (f_node[i][j].dp_flag == mini_dp_flag) {
                    if (f_node[i][j].score > max_score) {
                        max_score = f_node[i][j].score;
                        max_dis = f_node[i][j].dis_pen;
                        max_node = (line_node){i, j};
                        max_n = f_node[i][j].node_n;
                    } else if (f_node[i][j].score == max_score) {
                        if (f_node[i][j].dis_pen < max_dis) {
                            max_score = f_node[i][j].score;
                            max_dis = f_node[i][j].dis_pen;
                            max_node = (line_node){i, j};
                            max_n = f_node[i][j].node_n;
                        }
                    }
                }
            }
        }
    } else {
        f_node[right.x][right.y].from = head;//left;
        f_node[right.x][right.y].score = old_score;
        f_node[right.x][right.y].node_n = 1;
        frag_dp_update(f_node, a_msg, right.x, right.y, left.x+1, APP, AP, mini_dp_flag);
        max_score = f_node[right.x][right.y].score;
        max_node = f_node[right.x][right.y].from;
        max_n = f_node[right.x][right.y].node_n-1; // right is NOT inculding 
    }
    //backtrack
    line_node _right = max_node;
    int node_i = max_n-1;
    //while (_right.x != left.x)
    while (_right.x != head.x) {
        if (node_i < 0) { 
            fprintf(stderr, "\n[frag mini dp] %s node_i BUG 1.\n", APP->read_name);
            exit(0); 
        }
        line[node_i--] = _right;
        _right = f_node[_right.x][_right.y].from;
    }
    if (node_i >= 0) { fprintf(stderr, "\n[frag mini dp] %s node_i BUG 2.\n", APP->read_name); exit(0); }
    *de_score += (max_score-old_score);
    return max_n;
}

void frag_min_extend(frag_dp_node **f_node, aln_msg *a_msg,
                    int node_i, int aln_i,
                    int aln_min, int dp_flag,
                    lsat_aln_per_para *APP, lsat_aln_para *AP)
{
    int i, j, con_flag;
    int last_x = node_i, last_y = aln_i;
    // 0 -> node_i-1, node_i+1 -> node_n-1
    //from right to left 
    i = node_i - 1;
    while (i >= 0) {
        if (a_msg[i].n_aln > aln_min) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                get_fseed_dis(a_msg, i, j, node_i, aln_i, &con_flag, APP, AP);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH || con_flag == F_LONG_MISMATCH) {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        --i;
    }
    i = node_i + 1;
    while (i < APP->seed_out) {
        if (a_msg[i].n_aln > aln_min) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                get_fseed_dis(a_msg, node_i, aln_i, i, j, &con_flag, APP, AP);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH || con_flag == F_LONG_MISMATCH) {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        ++i;
    }
}

//Best Coverage and Connect
    //int frag_line_BCC(aln_msg *a_msg,
    //                  lsat_aln_per_para *APP, lsat_aln_para *AP,
    //                  line_node **line, int *line_end,
    //                  frag_dp_node ***f_node,
    //                  line_node **_line, int *_line_end)
    //{
    //    int i, j, k;
    //    int min_n = APP->min_thd, min_exist=0, min_num = 0;

    //    //dp init
    //    {
    //        for (i = 0; i < APP->seed_out; ++i)
    //        {
    //            if (a_msg[i].n_aln <= min_n)
    //            {
    //                frag_dp_init(*f_node, a_msg, i, (line_node){-1,0}, APP, AP, MIN_FLAG);
    //                min_exist = 1;
    //                ++min_num;
    //            }
    //            else
    //                frag_dp_init(*f_node, a_msg, i, (line_node){-1,0}, APP, AP, MULTI_FLAG);
    //        }
    //        //fraction XXX
    //        if (!min_exist || min_num * 10 < APP->seed_out)
    //        {
    //            for (i = 0; i < APP->seed_out; ++i)
    //            {
    //                for (j = 0; j < a_msg[i].n_aln; ++j) (*f_node)[i][j].dp_flag = MIN_FLAG;
    //            }
    //            min_n = APP->per_aln_n;
    //            min_exist = 1;
    //        }
    //    }
    //    //dp update 
    //    {
    //        //min extend, when min_n == PER_ALN_N: no need to extend
    //        if (min_n != APP->per_aln_n)
    //        {
    //            for (i = 0; i < APP->seed_out; ++i)
    //            {
    //                if (a_msg[i].n_aln <= min_n)
    //                {
    //                    for (j = 0; j < a_msg[i].n_aln; ++j)
    //                        frag_min_extend(*f_node, a_msg, i, j, min_n, MIN_FLAG, APP, AP);
    //                }
    //            }
    //        }
    //        //min update
    //        for (i = 1; i < APP->seed_out; ++i)
    //        {
    //            for (j = 0; j < a_msg[i].n_aln; ++j)
    //            {
    //                if ((*f_node)[i][j].dp_flag == MIN_FLAG)
    //                    frag_dp_update(*f_node, a_msg, i, j, 0/*update start pos*/, APP, AP, MIN_FLAG);
    //            }
    //        }
    //        //print
    //        /*printf("min-update:\n");
    //        for (i = 0; i < APP->seed_out; ++i)
    //        {
    //            for (j = 0; j < a_msg[i].n_aln; ++j)
    //                //if ((*f_node)[i][j].dp_flag == MIN_FLAG)
    //                    fprintf(stdout, "node:(%d %d)\t%d\t%d %d %lld\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].read_id, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, (long long)a_msg[i].at[j].offset, (*f_node)[i][j].from.x, (*f_node)[i][j].from.y, (*f_node)[i][j].score, FRAG_CON_STR[(*f_node)[i][j].match_flag], (*f_node)[i][j].dp_flag, (*f_node)[i][j].node_n);
    //        }*/
    //    }
    //    //backtrack
    //    {
    //        //find start node of backtrack
    //        int max_score=0, max_dis = 0, l_i = 0;
    //        line_node max_node = (line_node){-1,0};
    //        int node_i=0, mini_len, new_l=1;
    //        line_node last_n; int multi_l, max_multi = APP->seed_out > AP->res_mul_max ? AP->res_mul_max : APP->seed_out; //XXX
    //        line_node right, left;
    //        for (i = APP->seed_out-1; i >= 0; --i)
    //        {
    //            for (j = 0; j < a_msg[i].n_aln; ++j)
    //            {
    //                if ((*f_node)[i][j].backtrack_flag == DP_BACK_PEAK && (*f_node)[i][j].dp_flag == MIN_FLAG && (*f_node)[i][j].score > 0)
    //                {
    //                    if ((*f_node)[i][j].score > max_score)
    //                    {
    //                        max_score = (*f_node)[i][j].score;
    //                        max_node = (line_node){i,j};
    //                        max_dis = (*f_node)[i][j].dis_pen;
    //                    }
    //                    else if ((*f_node)[i][j].score == max_score)
    //                    {
    //                        if ((*f_node)[i][j].dis_pen < max_dis)
    //                        {
    //                            max_score = (*f_node)[i][j].score;
    //                            max_node = (line_node){i,j};
    //                            max_dis = (*f_node)[i][j].dis_pen;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //        //backtrack
    //        if (max_node.x == -1) return 0;

    //        int tri_x, tri_y;
    //        //ONLY one backtrack node.
    //        if (max_node.x < APP->seed_out-1)
    //        {
    //            //mini-dp with anchors
    //            mini_len = frag_mini_dp_line(*f_node, a_msg, APP, AP, max_node, (line_node){APP->seed_out, 0}, _line[0], 1, 0);

    //            if (mini_len > 0)
    //            {
    //                for (k = mini_len - 1; k >= 0; --k)
    //                    line[l_i][node_i++] = _line[0][k];	//_line[0][k] is the last node
    //            }

    //            line[l_i][node_i] = max_node;	//twice write, for the last multi-dp-line	//XXX
    //            //add trigger for multi-line
    //            last_n = (line_node){APP->seed_out, 0};
    //            for (k = mini_len; k >= 0; --k)
    //            {
    //                if (last_n.x - line[l_i][node_i-k].x > 2)
    //                {
    //                    tri_x = line[l_i][node_i-k].x;
    //                    tri_y = line[l_i][node_i-k].y;
    //                    //add multi-line
    //                    multi_l = frag_mini_dp_multi_line(*f_node, a_msg, APP, AP, line[l_i][node_i-k], last_n, _line, _line_end, max_multi);//, 0, 0);
    //                    for (i = 1; i < multi_l; ++i)
    //                    {
    //                        if (_line_end[i] > 0)
    //                        {
    //                            for (j = 0; j < _line_end[i]; ++j)
    //                                line[l_i+new_l][j] = _line[i][j];
    //                            line_end[l_i+new_l] = _line_end[i];
    //                            line[l_i+new_l][line_end[l_i+new_l]] = line[l_i][node_i-k];
    //                            line[l_i+new_l][line_end[l_i+new_l]+1] = last_n;
    //                            line[l_i+new_l][line_end[l_i+new_l]+2] = (line_node){0, 0};	//need trigger
    //                            if ((*f_node)[tri_x][tri_y].next_trigger_n == (*f_node)[tri_x][tri_y].trigger_m)
    //                            {
    //                                (*f_node)[tri_x][tri_y].trigger_m <<= 1;
    //                                (*f_node)[tri_x][tri_y].trigger = (int*)realloc((*f_node)[tri_x][tri_y].trigger,  (*f_node)[tri_x][tri_y].trigger_m * sizeof(int));
    //                            }
    //                            (*f_node)[tri_x][tri_y].trigger[(*f_node)[tri_x][tri_y].next_trigger_n++] = l_i+new_l;
    //                            ++new_l;
    //                        }
    //                    }
    //                }
    //                last_n = line[l_i][node_i-k];
    //            }
    //        }
    //        right = max_node;
    //        while (right.x != -1)
    //        {
    //            line[l_i][node_i++] = right;
    //            left = (*f_node)[right.x][right.y].from;

    //            if (left.x < right.x - 1 )//f_node[right.x][right.y].match_flag != F_MATCH)	//XXX if left.x < right.x-1, match_flag of right node couldn't be match?
    //            {
    //                //mini-dp with anchors
    //                mini_len = frag_mini_dp_line(*f_node, a_msg, APP, AP, left, right, _line[0], 1, 1);
    //                for (k = mini_len-1; k >= 0; --k)
    //                    line[l_i][node_i++] = _line[0][k];
    //                line[l_i][node_i] = left;	//twice write, for the last multi-dp-line	//XXX
    //                // add trigger multi-len
    //                last_n = right;
    //                for (k = mini_len; k >= 0; --k)
    //                {
    //                    if (last_n.x - line[l_i][node_i-k].x > 2)
    //                    {
    //                        //add multi-line
    //                        multi_l = frag_mini_dp_multi_line(*f_node, a_msg, APP, AP, line[l_i][node_i-k], last_n, _line, _line_end, max_multi);//, 0, 0);
    //                        for (i = 1; i < multi_l; ++i)
    //                        {
    //                            if (_line_end[i] > 0)
    //                            {
    //                                for (j = 0; j < _line_end[i]; ++j)
    //                                    line[l_i+new_l][j] = _line[i][j];
    //                                line_end[l_i+new_l] = _line_end[i];
    //                                line[l_i+new_l][line_end[l_i+new_l]] = line[l_i][node_i-k];
    //                                line[l_i+new_l][line_end[l_i+new_l]+1] = last_n;
    //                                line[l_i+new_l][line_end[l_i+new_l]+2] = (line_node){0, 0};	//need trigger
    //                                if ((*f_node)[last_n.x][last_n.y].next_trigger_n + (*f_node)[last_n.x][last_n.y].pre_trigger_n == (*f_node)[last_n.x][last_n.y].trigger_m)
    //                                {
    //                                    (*f_node)[last_n.x][last_n.y].trigger_m <<= 1;
    //                                    (*f_node)[last_n.x][last_n.y].trigger = (int*)realloc((*f_node)[last_n.x][last_n.y].trigger,  (*f_node)[last_n.x][last_n.y].trigger_m * sizeof(int));
    //                                }
    //                                (*f_node)[last_n.x][last_n.y].trigger[(*f_node)[last_n.x][last_n.y].next_trigger_n + (*f_node)[last_n.x][last_n.y].pre_trigger_n] = l_i+new_l;
    //                                ++(*f_node)[last_n.x][last_n.y].pre_trigger_n;
    //                                ++new_l;
    //                            }
    //                        }
    //                    }
    //                    last_n = line[l_i][node_i-k];
    //                }
    //            }
    //            right = left;
    //        }
    //        //invert line
    //        line_node tmp;
    //        for (k = 0; k < node_i/2; ++k)
    //        {
    //            tmp =line[l_i][k]; line[l_i][k] = line[l_i][node_i-k-1]; line[l_i][node_i-k-1] = tmp;
    //        }
    //        //print
    //        /*printf("whole-update:\n");
    //        for (i = 0; i < APP->seed_out; ++i)
    //        {
    //            for (j = 0; j < a_msg[i].n_aln; ++j)
    //            {
    //                fprintf(stdout, "node:(%d %d)\t%d\t%d %d %lld\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tBK-flag:%d\tnode_n:%d\n", i, j, a_msg[i].read_id, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, (long long)a_msg[i].at[j].offset, (*f_node)[i][j].from.x, (*f_node)[i][j].from.y, (*f_node)[i][j].score, FRAG_CON_STR[(*f_node)[i][j].match_flag], (*f_node)[i][j].dp_flag, (*f_node)[i][j].backtrack_flag, (*f_node)[i][j].node_n);
    //            }
    //        }*/
    //        line_end[l_i] = node_i;
    //        line[l_i][line_end[l_i]] = (line_node){-1,0};
    //        line[l_i][line_end[l_i]+1] = (line_node){APP->seed_out, 0};
    //        line[l_i][line_end[l_i]+2] = (line_node){1,1};	//trigger flag: (1,1) -> no needs for trigger
    //                                                        //				(0,0) -> need trigger
    //        
    //        /*for (i = 0; i < new_l; ++i)
    //        {
    //            fprintf(stdout, "line #%d:\n", i+1);
    //            for (j = 0; j < line_end[i]+3; ++j)
    //            {
    //                fprintf(stdout, "(%d,%d)\t", line[i][j].x, line[i][j].y);
    //            }
    //            fprintf(stdout, "\n");
    //        }*/
    //        return l_i+new_l;
    //    }
    //}

//new tree-generating and pruning
int frag_line_BCC(aln_msg *a_msg,
                  lsat_aln_per_para *APP, lsat_aln_para *AP,
                  line_node **line, int *line_end,
                  frag_dp_node ***f_node,
                  line_node **_line, int *_line_end, int line_n_max, int per_max_multi)
{
    int i, j, k;
    int min_n = APP->min_thd, min_exist=0, min_num = 0;

    { // dp init
        for (i = 0; i < APP->seed_out; ++i)
        {
            if (a_msg[i].n_aln <= min_n)
            {
                frag_dp_init(*f_node, a_msg, i, START_NODE, APP, AP, MIN_FLAG);
                min_exist = 1;
                ++min_num;
            }
            else
                frag_dp_init(*f_node, a_msg, i, START_NODE, APP, AP, MULTI_FLAG);
        }
        //fraction XXX
        if (!min_exist || min_num * 3 < APP->seed_out)
        {
            for (i = 0; i < APP->seed_out; ++i) {
                for (j = 0; j < a_msg[i].n_aln; ++j) (*f_node)[i][j].dp_flag = MIN_FLAG;
            }
            min_n = APP->per_aln_n;
            min_exist = 1;
        }
    }
    { // dp update
        //min extend, when min_n == PER_ALN_N: no need to extend
        if (min_n != APP->per_aln_n)
        {
            for (i = 0; i < APP->seed_out; ++i) {
                if (a_msg[i].n_aln <= min_n) {
                    for (j = 0; j < a_msg[i].n_aln; ++j)
                        frag_min_extend(*f_node, a_msg, i, j, min_n, MIN_FLAG, APP, AP);
                }
            }
        }
        //min update
        for (i = 1; i < APP->seed_out; ++i) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                if ((*f_node)[i][j].dp_flag == MIN_FLAG)
                    frag_dp_update(*f_node, a_msg, i, j, 0/*update start pos*/, APP, AP, MIN_FLAG);
            }
        }
    }
    { // backtrack
        //prune tree to find multi-nodes/roads for backtrack
        //XXX
        node_score *ns = node_init_score(line_n_max); //size of node_score
        line_node fa_node;
        for (i = APP->seed_out-1; i >=0; --i) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                if ((*f_node)[i][j].dp_flag == MIN_FLAG && (*f_node)[i][j].in_de == 0) //leaf node
                    branch_track_new(*f_node, i, j, ns); 
            }
        }
        int max_score, line_score; line_node max_node;
        int l_i = 0, new_l = 0, min_l=ns->node_n, o_l = min_l;
        int tri_x, tri_y, node_i, mini_len, multi_l;
        trig_node **inter_tri = (trig_node**)malloc(o_l* sizeof(trig_node*)); int *inter_n=(int*)malloc(o_l*sizeof(int));
        for (i = 0; i < o_l;++i)
            inter_tri[i] = (trig_node*)malloc(APP->seed_out * sizeof(trig_node));
        //multi-backtrack
        //extend for min-line AND set inter-trigger
        while (1) {
            //MAX score //XXX
            max_node = node_pop(ns, &line_score);    //这里得到的path,就是希望保留的
                                                     //总共ns->node_n条path,即min-line
            if (max_node.x == -1) break;
            //backtrack for max-node
            line_node last_n, right, left;
            //fill-up for right end
            node_i=0;
            inter_n[l_i]=0;
            if (max_node.x < APP->seed_out-1) {
                //XXX 是否有意义?
                mini_len = frag_mini_dp_line(*f_node, a_msg, APP, AP, max_node, (line_node){APP->seed_out, 0}, _line[0], &line_score, 1, 0);
                for (k = mini_len-1; k >= 0; --k) {
                    line[l_i][node_i++] = _line[0][k];
                    (*f_node)[_line[0][k].x][_line[0][k].y].dp_flag = MULTI_OUT; //不允许重复使用
                }
                line[l_i][node_i] = max_node;
                //set inter-trigger
                last_n = line[l_i][0];
                for (k = mini_len-1; k >= 0; --k) {
                    if (last_n.x - line[l_i][node_i-k].x > 2) {
                        // pre-trigger: (line[l_i][node_i-k].x , right.x)
                        (*f_node)[last_n.x][last_n.y].pre_trg = (line_node){line[l_i][node_i-k].x, last_n.x};
                        // 同一node处有多个trigger?XXX,即trg_n被重复赋值1或2
                        (*f_node)[last_n.x][last_n.y].trg_n |= 0x2;
                        inter_tri[l_i][inter_n[l_i]++] = (trig_node){(line_node){line[l_i][node_i-k].x, line[l_i][node_i-k].y}, (line_node){last_n.x, last_n.y}, 2};
                    }
                    last_n = line[l_i][node_i-k];
                }
            }
            right = max_node;
            // fill-up for every gap
            while (right.x != START_NODE.x) {
                line[l_i][node_i++] = right;
                left = (*f_node)[right.x][right.y].from;
                if (left.x < right.x - 1 )//f_node[right.x][right.y].match_flag != F_MATCH)	//XXX if left.x < right.x-1, match_flag of right node couldn't be match?
                {
                    //mini-dp with anchors
                    mini_len = frag_mini_dp_line(*f_node, a_msg, APP, AP, left, right, _line[0], &line_score, 1, 1);
                    for (k = mini_len-1; k >= 0; --k) {
                        line[l_i][node_i++] = _line[0][k];
                        (*f_node)[_line[0][k].x][_line[0][k].y].dp_flag = MULTI_OUT; //不允许重复使用
                    }
                    line[l_i][node_i] = left;	//twice write, for the last multi-dp-line	//XXX
                    //set inter-trigger
                    last_n = right;
                    for (k = mini_len; k >= 0; --k) {
                        if (last_n.x - line[l_i][node_i-k].x > 2) {
                            // pre-trigger: (line[l_i][node_i-k].x , right.x)
                            if (line[l_i][node_i-k].x == START_NODE.x) continue;
                            (*f_node)[last_n.x][last_n.y].pre_trg = (line_node){line[l_i][node_i-k].x, last_n.x};
                            (*f_node)[last_n.x][last_n.y].trg_n |= 0x2;
                            inter_tri[l_i][inter_n[l_i]++] = (trig_node){(line_node){line[l_i][node_i-k].x, line[l_i][node_i-k].y},(line_node){last_n.x, last_n.y}, 2};
                        }
                        last_n = line[l_i][node_i-k];
                    }
                }
                right = left;
            }
            //invert line
            line_node tmp;
            for (k = 0; k < node_i/2; ++k) {
                tmp =line[l_i][k]; line[l_i][k] = line[l_i][node_i-k-1]; line[l_i][node_i-k-1] = tmp;
            }
            line_end[l_i] = node_i;
            ///line[l_i][line_end[l_i]+2].x = line_score;
            L_LS(line, line_end, l_i) = line_score;
            l_i++;
        }
        /* line: [0 ~ end-1] -> node; 
                 [end] -> (left-b, right-b)/(start,end)
                 [end+1] -> (merged-flag, merge-head)/(start, end) // (-1,):NOT merged/ (-2,i):Merged and head is i/ (a,b): start from a to b
                 [end+2] -> (line-score, x) */
        // min-lines have been extended, and sorted by end-pos
        // set boundary and filter best/secondary
        line_set_bound(line, line_end, 0, &min_l, -1, APP->seed_out, inter_tri, inter_n, *f_node, a_msg, per_max_multi);
        for (i = 0; i < o_l; ++i) free(inter_tri[i]);
        free(inter_tri); free(inter_n);
        ///min_l = line_remove(line, line_end, 0, min_l);

        // set bound-trigger (-1, right)/ (left, APP->seed_out)
        for (i = 0; i < min_l; ++i) {
            if(L_MF(line, line_end, i) & L_INTER) continue;
            if (E_RB(line, line_end, i) == APP->seed_out && line[i][line_end[i]-1].x < APP->seed_out-2) {
                (*f_node)[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].next_trg = (line_node){line[i][line_end[i]-1].x, APP->seed_out};
                (*f_node)[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].trg_n |= 0x1;
            } 
            if (E_LB(line, line_end, i) == -1 && line[i][0].x > 1) {
                (*f_node)[line[i][0].x][line[i][0].y].pre_trg = (line_node){-1, line[i][0].x};
                (*f_node)[line[i][0].x][line[i][0].y].trg_n |= 0x2;
            }
        }
        //find path without endpoints between min-lines(merged-lines), similar to the trigger-line
        /*line_node left, right; int ll;
        right.x = APP->seed_out;
        for (k = 0; k < min_l; ++k) {
            if (L_MF(line, line_end, k) & L_INTER) continue;
            ///if (line[k][line_end[k]+1].x == -1) { //NO merged-line
            if (L_MF(line, line_end, k) & L_NMERG) {
                left.x = line[k][line_end[k]-1].x;
                // 两端的gap是 inter
                if (right.x != APP->seed_out && left.x != -1) {
                    multi_l = frag_mini_dp_multi_line(*f_node, a_msg, APP, AP, left.x, right.x, _line, _line_end, max_multi);//, 0, 0);
                    ll = new_l;
                    for (i = 0; i < multi_l; ++i) {
                        if (_line_end[i] > 0) {
                            for (j = 0; j < _line_end[i]; ++j) {
                                line[min_l+new_l][j] = _line[i][j];
                                (*f_node)[_line[i][j].x][_line[i][j].y].dp_flag = MULTI_OUT;
                            }
                            line_end[min_l+new_l] = _line_end[i];
                            //line[min_l+new_l][line_end[min_l+new_l]] = (line_node){left.x, right.x};
                            ++new_l;
                        }
                    }
                    line_set_bound1(line, line_end, min_l+ll, new_l-ll, left.x, right.x);
                    new_l = ll+line_remove(line, line_end, min_l+ll, new_l-ll);
                }
                right.x = line[k][0].x;
            ///} else if (line[k][line_end[k]+1].x != -2) { //merged and head
            } else if (L_MF(line, line_end, k) & L_MERGH) {
                ///left.x = line[k][line_end[k]+1].y;
                left.x = L_RB(line, line_end, k);
                if (right.x != APP->seed_out && left.x != -1) {
                    multi_l = frag_mini_dp_multi_line(*f_node, a_msg, APP, AP, left.x, right.x, _line, _line_end, max_multi);//, 0, 0);
                    ll = new_l;
                    for (i = 0; i < multi_l; ++i) {
                        if (_line_end[i] > 0) {
                            for (j = 0; j < _line_end[i]; ++j) {
                                line[min_l+new_l][j] = _line[i][j];
                                (*f_node)[_line[i][j].x][_line[i][j].y].dp_flag = MULTI_OUT;
                            }
                            line_end[min_l+new_l] = _line_end[i];
                            ++new_l;
                        }
                    }
                    line_set_bound1(line, line_end, min_l+ll, new_l-ll, left.x, right.x);
                    new_l = ll + line_remove(line, line_end, min_l+ll, new_l-ll);
                }
                ///right.x = line[k][line_end[k]+1].x;
                right.x = L_RB(line, line_end, k);
            } // merged and body: continue
        }*/
        node_free_score(ns);
        return min_l+new_l;
    }
}

int frag_line_remain(aln_reg *re_reg)
{
}

//_dp
    //all valid aln 
    //@para:
    //	a_msg:	struct of seeds' aln result
    //	n_seed:	whole number of seeds that have at least 1 aln-res, and less than 100
    //	f_msg:	be used to store frag-msg
    //	line_n:	number of lines
    //	line_m:	max number of lines allowed in mem
    //int _frag_dp_path(aln_msg *a_msg, 
    //        lsat_aln_per_para *APP, lsat_aln_para *AP,
    //        frag_msg **f_msg, 
    //        int *line_n/*return, line_num*/, int *line_tri, int *line_m/*return, line mem size*/, 
    //        line_node **line, int *line_end, 
    //        frag_dp_node ***f_node, 
    //        line_node **_line, int *_line_end,
    //        int aln_type)
    //{
    //    int i, l;
    //    //DP
    //    int l_n;
    //    if (aln_type == 1) // best coverage and connect
    //        l_n = frag_line_BCC(a_msg, APP, AP, line, line_end, f_node, _line, _line_end, 10);// for best coverage and connect case: l_n equals 1.
    //    else l_n = 0; // all valid
    //    if (l_n == 0) return 0;
    //
    //    /*for (i = 0; i < l_n; ++i)
    //      {
    //      printf("%d:\t", i+1);
    //      for (l = 0; l < line_end[i]; ++l)
    //      printf("(%d, %d)\t", a_msg[line[i][l].x].read_id, line[i][l].y);
    //      printf("\n");
    //      }*/
    //
    //    int frag_num;
    //    int cur_x , cur_y , pre_x, pre_y;
    //
    //    //for new frags
    //    //  cur_line: cur line index, cur_num:  whole lines, *line_num: mem of *f_msg
    //    int cur_line=0;//, cur_num=1, _m_len;
    //    int left_bound, right_bound;
    //
    //    // calcu whole number of line
    //    (*line_n) = 1;
    //    for (i = line_end[0]-1; i > 0; --i)
    //    {
    //        if ((*f_node)[line[0][i].x][line[0][i].y].match_flag == F_CHR_DIF || (*f_node)[line[0][i].x][line[0][i].y].match_flag == F_UNCONNECT)
    //            ++(*line_n);
    //    }
    //    if (*line_n > (*line_m))
    //    {
    //        (*f_msg) = (frag_msg*)realloc(*f_msg, *line_n * sizeof(frag_msg));
    //        for (i = (*line_m); i < *line_n; ++i)
    //        {
    //            frag_init_msg((*f_msg)+i, (*f_msg)->frag_max); 
    //            frag_copy_msg(*f_msg, (*f_msg)+i);
    //        }
    //        if ((*f_msg) == NULL) { fprintf(stderr, "\n[frag_dp_path] Not enough memory.(line_m: %d)\n", *line_n); exit(0); }
    //        (*line_m)= *line_n;
    //        fprintf(stderr, "line-num: %d\t", *line_n);
    //    }
    //
    //    for (l = 0; l < l_n; ++l)
    //    {
    //        frag_num=0;
    //        pre_x = line[l][line_end[l]-1].x;
    //        pre_y = line[l][line_end[l]-1].y;
    //        //fprintf(stdout, "left: %d, right: %d\n", line[l][line_end[l]].x, line[l][line_end[l]+1].x);
    //        //right bound
    //        right_bound = ((line[l][line_end[l]+1].x == APP->seed_out) ? (APP->seed_all+1) : (a_msg[line[l][line_end[l]+1].x].read_id));
    //        //first end
    //        frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num);
    //        ((*f_msg)+cur_line)->frag_right_bound = right_bound;
    //        line_tri[cur_line] = 1;
    //        for (i = line_end[l]-1; i > 0; --i)
    //        {
    //            cur_x = pre_x; cur_y = pre_y;
    //            pre_x = line[l][i-1].x; pre_y = line[l][i-1].y;
    //            //INS
    //            if ((*f_node)[cur_x][cur_y].match_flag == F_INSERT)
    //            {
    //                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num);
    //                ++frag_num;
    //                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num);
    //            }
    //            //DEL
    //            else if ((*f_node)[cur_x][cur_y].match_flag == F_DELETE)
    //            {
    //                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num);
    //                ++frag_num;
    //                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num);
    //            }
    //            //MIS	XXX
    //            else if ((*f_node)[cur_x][cur_y].match_flag == F_MISMATCH || (*f_node)[cur_x][cur_y].match_flag == F_LONG_MISMATCH)
    //            {
    //                //XXX take mis-match as NEW-FLAG case
    //                //frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
    //                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num);
    //                ++frag_num;
    //                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num);
    //            }
    //            else if ((*f_node)[cur_x][cur_y].match_flag == F_MATCH)	//: MATCH
    //            {
    //                frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num);
    //            }
    //            else if ((*f_node)[cur_x][cur_y].match_flag == F_CHR_DIF) // : INV
    //            {
    //                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num);
    //                left_bound = a_msg[pre_x].read_id;
    //                ((*f_msg)+cur_line)->frag_left_bound = left_bound;
    //
    //                right_bound = a_msg[cur_x].read_id;
    //                ++cur_line;
    //                if (cur_line >= (*line_m)) fprintf(stderr, "[frag dp path] line number error.\n");
    //                frag_num = 0;
    //                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num);
    //                ((*f_msg)+cur_line)->frag_right_bound = right_bound;
    //                line_tri[cur_line] = 1;
    //            }
    //            else if ((*f_node)[cur_x][cur_y].match_flag == F_UNCONNECT) // : Alu ...
    //            {
    //                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num);
    //                left_bound = a_msg[pre_x].read_id;
    //                ((*f_msg)+cur_line)->frag_left_bound = left_bound;
    //
    //                right_bound = a_msg[cur_x].read_id;
    //                ++cur_line;
    //                if (cur_line >= (*line_m)) fprintf(stderr, "[frag dp path] line number error.\n");
    //                frag_num = 0;
    //                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num);
    //                ((*f_msg)+cur_line)->frag_right_bound = right_bound;
    //                line_tri[cur_line] = 1;
    //            }
    //            else {fprintf(stderr, "[frag dp path] Error: Unknown flag, \"%d\"\n", (*f_node)[cur_x][cur_y].match_flag); exit(0);}
    //        }
    //        //last start
    //        cur_x = line[l][0].x; cur_y = line[l][0].y;
    //        frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num);
    //        left_bound = ((line[l][line_end[l]].x == -1) ? 0 : (a_msg[line[l][line_end[l]].x].read_id));
    //        ((*f_msg)+cur_line)->frag_left_bound = left_bound;
    //    }
    //
    //    /*int k, s_i, a_i;
    //      fprintf(stdout, "%d line(s)\n", *line_n);
    //      for (i = 0; i < *line_n; ++i)
    //      {
    //      fprintf(stdout, "line#: %d:\tleft: %d  right: %d\n", i+1, ((*f_msg)+i)->frag_left_bound, ((*f_msg)+i)->frag_right_bound);
    //      for (j = 0; j < ((*f_msg)+i)->frag_num; ++j)
    //      {
    //      fprintf(stdout, "\tfrag#: %d\n", j+1);
    //      for (k = 0; k < ((*f_msg)+i)->fa_msg[j].seed_num; ++k)
    //      {
    //      s_i = ((*f_msg)+i)->fa_msg[j].seed_i[k];
    //      a_i = ((*f_msg)+i)->fa_msg[j].seed_aln_i[k];
    //      fprintf(stdout, "\t\t%d %d %d %lld\n", a_msg[s_i].read_id, a_msg[s_i].at[a_i].nsrand, a_msg[s_i].at[a_i].chr, (long long)a_msg[s_i].at[a_i].offset);
    //      }
    //      }
    //      }*/
    //    return 1;
    //}
int frag_dp_path(aln_msg *a_msg, frag_msg **f_msg,
                 lsat_aln_per_para *APP, lsat_aln_para *AP,
                 int *l_n, int *line_m,
                 line_node **line, int *line_end,
                 frag_dp_node ***f_node,
                 line_node **_line, int *_line_end)
{
    int line_n = *l_n;
    if (line_n == 0) return 0;
	int i, l;
    
#ifdef __DEBUG__
    for (i = 0; i < line_n; ++i) {
        printf("%d(%d,%d score: %d):\t", i+1, L_MF(line, line_end, i), L_MH(line, line_end,i), L_LS(line, line_end, i));
        for (l = 0; l < line_end[i]; ++l)
            printf("(%d, %d)\t", a_msg[line[i][l].x].read_id-1, line[i][l].y);
            //printf("(%d, %d)\t", line[i][l].x, line[i][l].y);
        printf("\n");
    }
#endif

	int frag_num;
	int cur_x , cur_y , pre_x, pre_y;

	int left_bound, right_bound;

    if (line_n > (*line_m)) {
        (*f_msg) = (frag_msg*)realloc(*f_msg, line_n * sizeof(frag_msg));
        for (i = (*line_m); i < line_n; ++i)
        {
            frag_init_msg((*f_msg)+i, (*f_msg)->frag_max); 
            frag_copy_msg(*f_msg, (*f_msg)+i);
        }
        if ((*f_msg) == NULL) { fprintf(stderr, "\n[frag_dp_path] Not enough memory.(line_m: %d)\n", line_n); exit(0); }
        (*line_m)= line_n;
        fprintf(stderr, "line-num: %d\t", line_n);
    }

    int j=0;
	for (l = 0; l < line_n; ++l) {
        //set merged msg
        ((*f_msg)+j)->merg_msg.x = L_MF(line,line_end,l);
        ((*f_msg)+j)->merg_msg.y = L_MH(line,line_end,l);
            
        frag_num=0;
        pre_x = line[l][line_end[l]-1].x; pre_y = line[l][line_end[l]-1].y;
        //fprintf(stdout, "left: %d, right: %d\n", line[l][line_end[l]].x, line[l][line_end[l]+1].x);
		//right bound
        ///right_bound = ((line[l][line_end[l]].y == APP->seed_out) ? (APP->seed_all+1) : (a_msg[line[l][line_end[l]].y].read_id));
        right_bound = ((E_RB(line, line_end, l) == APP->seed_out) ? (APP->seed_all+1) : a_msg[E_RB(line, line_end, l)].read_id);
        //MIS-MATCH
        //first end
		frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num); //XXX set merged msg
		if ((*f_node)[pre_x][pre_y].trg_n)	//set trigger
            frag_trg_set((*f_node)[pre_x][pre_y], (*f_msg)+j, frag_num);
        ((*f_msg)+j)->frag_right_bound = right_bound;
		for (i = line_end[l]-1; i > 0; --i) {
			cur_x = pre_x; cur_y = pre_y;
			pre_x = line[l][i-1].x; pre_y = line[l][i-1].y;
			//INS
			if ((*f_node)[cur_x][cur_y].match_flag == F_INSERT)
			{
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
				if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
                    frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num);
            }
			//DEL
			else if ((*f_node)[cur_x][cur_y].match_flag == F_DELETE) {
				frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
				if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
                    frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num);
			}
			//MIS	XXX
			else if ((*f_node)[cur_x][cur_y].match_flag == F_MISMATCH || (*f_node)[cur_x][cur_y].match_flag == F_LONG_MISMATCH) {
				//XXX take mis-match as NEW-FLAG case
				//frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
				frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
				if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
                    frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num);
				//XXX find new flag in mis-match seeds
				//dp-line for INV/TRS
            } else if ((*f_node)[cur_x][cur_y].match_flag == F_MATCH)	//: MATCH
			{
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+j, frag_num);
			} else { fprintf(stderr, "\n[frag dp path] Error: Unknown flag, \"%d\"\n", (*f_node)[cur_x][cur_y].match_flag); exit(0);}
		}
		//last start
		cur_x = line[l][0].x; cur_y = line[l][0].y;
		frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
        if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
            frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
        ///left_bound = ((line[l][line_end[l]].x == -1) ? 0 : (a_msg[line[l][line_end[l]].x].read_id));
        left_bound = ((E_LB(line, line_end, l) == -1) ? 0 : a_msg[E_LB(line, line_end, l)].read_id);
		//MIS-MATCH
        ((*f_msg)+j)->frag_left_bound = left_bound;
#ifdef __DEBUG__
        fprintf(stdout, "#%d: left %d, right %d\n", j+1, left_bound, right_bound);
#endif
        ++j;
	}
    *l_n = j;
	
	/*int k, s_i, a_i;
	fprintf(stdout, "%d line(s)\n", *line_n);
	for (i = 0; i < *line_n; ++i)
	{
		fprintf(stdout, "line#: %d:\tleft: %d  right: %d\n", i+1, ((*f_msg)+i)->frag_left_bound, ((*f_msg)+i)->frag_right_bound);
		for (j = 0; j < ((*f_msg)+i)->frag_num; ++j)
		{
			fprintf(stdout, "\tfrag#: %d\n", j+1);
			for (k = 0; k < ((*f_msg)+i)->fa_msg[j].seed_num; ++k)
			{
				s_i = ((*f_msg)+i)->fa_msg[j].seed_i[k];
				a_i = ((*f_msg)+i)->fa_msg[j].seed_aln_i[k];
				fprintf(stdout, "\t\t%d %d %d %lld\n", a_msg[s_i].read_id, a_msg[s_i].at[a_i].nsrand, a_msg[s_i].at[a_i].chr, (long long)a_msg[s_i].at[a_i].offset);
			}
		}
	}*/
    return 1;
}

void lsat_unmap(char *read_name)
{
    fprintf(stdout, "%s\t*\t*\t*\t*\n", read_name);
}

int frag_sam_cluster(const char *read_prefix, char *seed_result, seed_msg *s_msg, bwt_t *bwt, bntseq_t *bns, uint8_t *pac, lsat_aln_per_para *APP, lsat_aln_para *AP)
{
    int read_n/*1-based*/, seed_out, i, j;
    int seed_id;

    sam_msg *m_msg;
    samfile_t *samf = 0;

    aln_msg *a_msg; 
    frag_msg *f_msg; int line_n=0, tri_n, line_m=1/*line size*/;
    frag_dp_node ***f_node;
    uint32_t *hash_num;
    uint64_t **hash_node;

    gzFile readfp; kseq_t *read_seq_t; char *read_seq;

    //alloc mem and initialization
        m_msg = sam_init_msg();
        a_msg = aln_init_msg(s_msg->seed_max, AP->per_aln_m);
        f_node = fnode_alloc(s_msg->seed_max+2, AP->per_aln_m);

        int line_n_max = s_msg->seed_max * AP->per_aln_m;
        int line_len_max = s_msg->seed_max;
        line_node **line = (line_node**)malloc(line_n_max * sizeof(line_node*));
        int *line_end = (int*)malloc(line_n_max * sizeof(int));
        line_node **_line = (line_node**)malloc(line_n_max * sizeof(line_node*));
        int *_line_end = (int*)malloc(line_n_max * sizeof(int));
        for (i = 0; i < line_n_max; ++i)
        {
            line[i] = (line_node*)malloc((line_len_max+L_EXTRA) * sizeof(line_node));
            _line[i] = (line_node*)malloc((line_len_max+L_EXTRA) * sizeof(line_node));
        }

        if (line == NULL || _line == NULL) { fprintf(stderr, "\n[frag_dp_path] Not enougy memory.\n"); exit(0); }

        //XXX 
        f_msg = (frag_msg*)malloc(sizeof(frag_msg));	
        frag_init_msg(&f_msg[0], s_msg->seed_max);


        //alloc mem for hash mapping
        int key_len = 2;
        int hash_size = (int)pow(NT_N, key_len);
        hash_num = (uint32_t*)calloc(hash_size, sizeof(uint32_t));	//16 = pow(4, 2)
        hash_node = (uint64_t**)calloc(hash_size, sizeof(uint64_t*));

    if ((samf = samopen(seed_result, "r")) == 0) {
        fprintf(stderr, "\n[lsat_aln] Can't open seed result sam file %s.\n", seed_result);
        exit(-1);
    }
    if (samf->header == 0) {
        fprintf(stderr, "\n[lsat_aln] Can't read the header of result sam file %s.\n", seed_result);
        exit(-1);
    }
    int r;
    seed_id = seed_out = 0;
    read_n = 1;

    readfp = gzopen(read_prefix, "r");
    read_seq_t = kseq_init(readfp);
    while (read_n <= s_msg->read_all) { // get seeds' aln-msg of every read
        init_aln_per_para(APP, s_msg, read_n);
        if ((r = sam_read1(samf->x.tamr, samf->header, m_msg, APP->per_aln_n)) < 0) break;
        if (APP->per_aln_n > AP->per_aln_m)
        {
            aln_realc_msg(a_msg, s_msg->seed_max, AP->per_aln_m, APP->per_aln_n);
            fnode_realc(f_node, s_msg->seed_max+2, AP->per_aln_m, APP->per_aln_n);
            AP->per_aln_m = APP->per_aln_n;
        }
        ++seed_id;
        if (r > 0) {
            ++seed_out;
            for (i = 0; i < m_msg->sam_n; ++i) 
                setAmsg(a_msg, seed_out, i+1, seed_id - s_msg->seed_all[read_n-1], m_msg->sam[i]);
        }
        if (seed_id == s_msg->seed_all[read_n]) { // get a whole-read
            if (kseq_read(read_seq_t) < 0) { fprintf(stderr, "\n[lsat_aln] Read file ERROR.\n"); exit(-1); }
            // XXX no seed is 'setAmsg'ed
            read_seq = read_seq_t->seq.s;

            APP->seed_out = seed_out;
            //fnode_init(*f_node, s_msg->seed_max+2, AP->per_aln_m);
            
            //int per_max_multi = APP->seed_out > AP->res_mul_max ? AP->res_mul_max : APP->seed_out; // XXX
			int per_max_multi = AP->res_mul_max;
            aln_res *a_res = aln_init_res(1, APP->read_len); // line_n
            aln_reg *a_reg = aln_init_reg(APP->read_len);
            // main line process
            strcpy(READ_NAME, APP->read_name);
            line_n = frag_line_BCC(a_msg, APP, AP, line, line_end, f_node, _line, _line_end, line_n_max, per_max_multi);
            if (line_n == 0) 
                lsat_unmap(s_msg->read_name[read_n]);
            else {
                frag_dp_path(a_msg, &f_msg, APP, AP, &line_n, &line_m, line, line_end, f_node, _line, _line_end);
                frag_check(a_msg, &f_msg, a_res, bns, pac, read_prefix, read_seq, APP, AP, line_n, &hash_num, &hash_node);
                aln_res_output(a_res, a_reg, APP);
                // trigger-line
                for (i=0; i<line_n; ++i) {
                    if (a_res->la[i].merg_msg.x == 1) {
                        for (j=0; j<a_res->la[i].trg_n; ++j) {
                            tri_n = trg_dp_line(*f_node, a_msg, &f_msg, APP, AP, a_res->la[i].trg[j].x, a_res->la[i].trg[j].y, line, line_end, _line, _line_end, per_max_multi);
                            if (tri_n == 0) continue;
                            frag_dp_path(a_msg, &f_msg, APP, AP, &tri_n, &line_m, line, line_end, f_node, _line, _line_end);
                            aln_res *tri_res = aln_init_res(1, APP->read_len);
                            frag_check(a_msg, &f_msg, tri_res, bns, pac, read_prefix, read_seq, APP, AP, tri_n, &hash_num, &hash_node);
                            aln_res_output(tri_res, a_reg, APP); aln_res_free(tri_res);
                        }
                    }
                }
                aln_res_free(a_res); aln_free_reg(a_reg);
            }
            seed_out = 0;
            ++read_n;
        }
    }
    //free variables and close file handles
    for (i = 0; i < m_msg->sam_m; ++i) {
        if (m_msg->sam[i].cigar_s->s) free(m_msg->sam[i].cigar_s->s);
        free(m_msg->sam[i].cigar_s);
    } 
    free(m_msg->sam); free(m_msg);
    aln_free_msg(a_msg, s_msg->seed_max, AP->per_aln_m);
    fnode_free(f_node, s_msg->seed_max+2, AP->per_aln_m);
    
    for (i = 0; i < line_n_max; ++i) {
        free(line[i]); free(_line[i]);
    }
    free(line); free(_line); free(line_end); /*free(line_tri);*/ free(_line_end);
    //hash map
    for (i = 0; i < hash_size; ++i) free(hash_node[i]);
    free(hash_node); free(hash_num); 

    frag_free_msg(f_msg, line_m);
    gzclose(readfp); kseq_destroy(read_seq_t);
    samclose(samf);
    return 0;
}

// for GEM
int frag_map_cluster(const char *read_prefix, char *seed_result, seed_msg *s_msg, bwt_t *bwt, bntseq_t *bns, uint8_t *pac, lsat_aln_per_para *APP, lsat_aln_para *AP)
{
    int read_n/*1-based*/, seed_out, i, j;
    int seed_id;

    map_msg *m_msg;
    FILE *mapf;

    aln_msg *a_msg; 
    frag_msg *f_msg; int line_n=0, tri_n, line_m=1/*line size*/;
    frag_dp_node ***f_node;
    uint32_t *hash_num;
    uint64_t **hash_node;

    gzFile readfp; kseq_t *read_seq_t; char *read_seq;

    //alloc mem and initialization
        m_msg = map_init_msg();
        a_msg = aln_init_msg(s_msg->seed_max, AP->per_aln_m);
        f_node = fnode_alloc(s_msg->seed_max+2, AP->per_aln_m);

        int line_n_max = s_msg->seed_max * AP->per_aln_m;
        int line_len_max = s_msg->seed_max;
        line_node **line = (line_node**)malloc(line_n_max * sizeof(line_node*));
        int *line_end = (int*)malloc(line_n_max * sizeof(int));
        line_node **_line = (line_node**)malloc(line_n_max * sizeof(line_node*));
        int *_line_end = (int*)malloc(line_n_max * sizeof(int));
        for (i = 0; i < line_n_max; ++i) {
            line[i] = (line_node*)malloc((line_len_max+L_EXTRA) * sizeof(line_node));
            _line[i] = (line_node*)malloc((line_len_max+L_EXTRA) * sizeof(line_node));
        }

        if (line == NULL || _line == NULL) { fprintf(stderr, "\n[frag_dp_path] Not enougy memory.\n"); exit(0); }

        //XXX 
        f_msg = (frag_msg*)malloc(sizeof(frag_msg));	
        frag_init_msg(&f_msg[0], s_msg->seed_max);


        //alloc mem for hash mapping
        int key_len = 2;
        int hash_size = (int)pow(NT_N, key_len);
        hash_num = (uint32_t*)calloc(hash_size, sizeof(uint32_t));	//16 = pow(4, 2)
        hash_node = (uint64_t**)calloc(hash_size, sizeof(uint64_t*));

    if ((mapf = fopen(seed_result, "r")) == NULL) {
        fprintf(stderr, "\n[lsat_aln] Can't open seed result 'map' file %s.\n", seed_result);
        exit(-1);
    }

    int r;
    seed_id = seed_out = 0;
    read_n = 1;

    readfp = gzopen(read_prefix, "r");
    read_seq_t = kseq_init(readfp);
    while (read_n <= s_msg->read_all) { // get seeds' aln-msg of every read
        init_aln_per_para(APP, s_msg, read_n);
        if ((r = gem_map_read(mapf, m_msg, APP->per_aln_n)) < 0) break;
        if (APP->per_aln_n > AP->per_aln_m) {
            aln_realc_msg(a_msg, s_msg->seed_max, AP->per_aln_m, APP->per_aln_n);
            fnode_realc(f_node, s_msg->seed_max+2, AP->per_aln_m, APP->per_aln_n);
            AP->per_aln_m = APP->per_aln_n;
        }
        ++seed_id;
        if (r > 0) {
            ++seed_out;
            for (i = 0; i < m_msg->map_n; ++i) 
                set_aln_msg(a_msg, seed_out, i+1, seed_id - s_msg->seed_all[read_n-1], m_msg->map[i]);
        }
        if (seed_id == s_msg->seed_all[read_n]) { // get a whole-read
            if (kseq_read(read_seq_t) < 0) { fprintf(stderr, "\n[lsat_aln] Read file ERROR.\n"); exit(-1); }
            // XXX no seed is 'setAmsg'ed
            read_seq = read_seq_t->seq.s;

            APP->seed_out = seed_out;
            //fnode_init(*f_node, s_msg->seed_max+2, AP->per_aln_m);
            
            //int per_max_multi = APP->seed_out > AP->res_mul_max ? AP->res_mul_max : APP->seed_out; // XXX
			int per_max_multi = AP->res_mul_max;
            aln_res *a_res = aln_init_res(1, APP->read_len);
			aln_res *remain_res = aln_init_res(1, APP->read_len);
            aln_reg *a_reg = aln_init_reg(APP->read_len);
            // main line process
            strcpy(READ_NAME, APP->read_name);
            line_n = frag_line_BCC(a_msg, APP, AP, line, line_end, f_node, _line, _line_end, line_n_max, per_max_multi);
            if (line_n == 0) 
                lsat_unmap(s_msg->read_name[read_n]);
            else {
                frag_dp_path(a_msg, &f_msg, APP, AP, &line_n, &line_m, line, line_end, f_node, _line, _line_end);
                frag_check(a_msg, &f_msg, a_res, bns, pac, read_prefix, read_seq, APP, AP, line_n, &hash_num, &hash_node);
                aln_res_output(a_res, a_reg, APP);
				/*
                // trigger-line
                for (i=0; i<line_n; ++i) {
                    if (a_res->la[i].merg_msg.x == 1) {
                        for (j=0; j<a_res->la[i].trg_n; ++j) {
                            tri_n = trg_dp_line(*f_node, a_msg, &f_msg, APP, AP, a_res->la[i].trg[j].x, a_res->la[i].trg[j].y, line, line_end, _line, _line_end, per_max_multi);
                            if (tri_n == 0) continue;
                            frag_dp_path(a_msg, &f_msg, APP, AP, &tri_n, &line_m, line, line_end, f_node, _line, _line_end);
                            aln_res *tri_res = aln_init_res(1, APP->read_len);
                            frag_check(a_msg, &f_msg, tri_res, bns, pac, read_prefix, read_seq, APP, AP, tri_n, &hash_num, &hash_node);
                            aln_res_output(tri_res, a_reg, APP); aln_res_free(tri_res);
                        }
                    }
                }*/
				// re-DP for remian regions with long-seed-reuslts
				frag_line_remain(a_reg);
                // bwt aln for remain regions 
                bwt_aln_remain(a_reg, remain_res, bwt, bns, pac, read_seq, APP, AP);
                aln_res_output(remain_res, 0, APP);
            }
			// free 
			aln_res_free(a_res); aln_res_free(remain_res);
			aln_free_reg(a_reg);
            seed_out = 0;
            //if (read_n % 1000 == 0) 
			fprintf(stderr, "%16d reads have been aligned.\n", read_n);
            ++read_n;
        }
    }
    //free variables and close file handles
    map_free_msg(m_msg);
    aln_free_msg(a_msg, s_msg->seed_max, AP->per_aln_m);
    fnode_free(f_node, s_msg->seed_max+2, AP->per_aln_m);
    
    for (i = 0; i < line_n_max; ++i) {
        free(line[i]); free(_line[i]);
    }
    free(line); free(_line); free(line_end); /*free(line_tri);*/ free(_line_end);
    //hash map
    for (i = 0; i < hash_size; ++i) free(hash_node[i]);
    free(hash_node); free(hash_num); 

    frag_free_msg(f_msg, line_m);
    gzclose(readfp); kseq_destroy(read_seq_t);
    fclose(mapf);
    return 0;
}

int lsat_gem(const char *ref_prefix, const char *read_prefix, lsat_aln_per_para *APP)
{
    char cmd[1024];
    sprintf(cmd, "./gem_map.sh %s %s.seed -d %d", ref_prefix, read_prefix, lsat_per_aln[0]); // lsat_per_aln[0]: max number of loctions per seed
    fprintf(stderr, "[lsat_aln] Executing gem-mapper ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lsat_aln] Seeding undone, gem-mapper exit abnormally.\n"); exit(0); }
    fprintf(stderr, "done.\n");
    return 0;
}


int lsat_bwa(const char *ref_prefix, const char *read_prefix)
{
    char cmd[1024];
    sprintf(cmd, "./bwa_aln.sh %s %s.seed", ref_prefix, read_prefix);
    fprintf(stderr, "[lsat_aln] Executing bwa aln ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lsat_aln] Seeding undone, bwa aln exit abnormally.\n"); exit(0); }
    fprintf(stderr, "done.\n");
    return 0;
}

int lsat_soap2dp(const char *ref_prefix, const char *read_prefix)
{
    char cmd[1024];
    sprintf(cmd, "./soap2dp_aln.sh %s %s", ref_prefix, read_prefix);
    fprintf(stderr, "[lsat_aln] Executing soap2dp aln ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lsat_aln] Seeding undone, soap2dp aln exit abnormalloy.\n"); exit(0); }
    fprintf(stderr, "done.\n");
    return 0;
}

int lsat_aln_core(const char *ref_prefix, const char *read_prefix, int seed_info, int seed_program, int no_seed_aln, char *seed_result, lsat_aln_per_para *APP, lsat_aln_para *AP)
{
    clock_t t = clock();
    seed_msg *s_msg;

    /* split-seeding */
    s_msg = seed_init_msg();
    if (seed_info) split_seed_info(read_prefix, s_msg);
    else split_seed(read_prefix, s_msg);

    if (!strcmp(seed_result, ""))
    {
        strcpy(seed_result, read_prefix);
        if (seed_program == 0)
            strcat(seed_result, ".seed.gem.map");
        else if (seed_program == 1)
            strcat(seed_result, ".seed.bwa.sam");
        else if (seed_program == 2)
            strcat(seed_result, ".seed.out.0");
        else { fprintf(stderr, "[lsat_aln] Unknown seeding program option.\n"); return lsat_aln_usage(); }
    }
    //excute soap2-dp program
    if (!no_seed_aln) 
    {
        if (seed_program == 0) lsat_gem(ref_prefix, read_prefix, APP);
        else if (seed_program == 1) lsat_bwa(ref_prefix, read_prefix);
        else if (seed_program == 2) lsat_soap2dp(ref_prefix, read_prefix);
        else { fprintf(stderr, "[lsat_aln] Unknown seeding program option.\n"); return lsat_aln_usage(); }
    }

    /* frag-clustering */
    /* load index */
    fprintf(stderr, "[lsat_aln] Restoring ref-indices ... ");
    bwt_t *bwt; bntseq_t *bns;
    char *str = (char*)calloc(strlen(ref_prefix)+5, 1);
    strcpy(str, ref_prefix); strcat(str, ".bwt"); bwt = bwt_restore_bwt(str);
    strcpy(str, ref_prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
    free(str);

    bns = bns_restore(ref_prefix);
    uint8_t *pac = (uint8_t*)calloc(bns->l_pac/4+1, 1);
    fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);	fprintf(stderr, "done.\n");

    fprintf(stderr, "[lsat_aln] Clustering frag ... ");
    if (seed_program == 0)
        frag_map_cluster(read_prefix, seed_result, s_msg, bwt,  bns, pac, APP, AP);	
    else
        frag_sam_cluster(read_prefix, seed_result, s_msg, bwt, bns, pac, APP, AP);	
    fprintf(stderr, "done.\n");
    fprintf(stderr, "[lsat_aln] Time Consupmtion: %.3f sec.\n", (float)(clock() -t )/CLOCKS_PER_SEC);

    seed_free_msg(s_msg);
    bwt_destroy(bwt); bns_destroy(bns); free(pac);
    return 0;
}

int lsat_aln(int argc, char *argv[])
{
    int c;
    char *ref, *read;
    // parameters
    lsat_aln_per_para *APP = (lsat_aln_per_para*)calloc(1, sizeof(lsat_aln_per_para));
    lsat_aln_para *AP = (lsat_aln_para*)calloc(1, sizeof(lsat_aln_para));
    init_aln_para(AP);
    int no_seed_aln=0, seed_info=0, aln_type=1, seed_program=0;
	int i, seed_len=0, seed_inv=0, per_aln=0, min_thd=0;
    char seed_result_f[1024]="", aln_result_f[1024]="";

    while ((c =getopt(argc, argv, "t:m:M:O:E:S:r:V:g:s:l:v:p:d:o:NIA:")) >= 0)
    {
        switch (c)
        {
			case 't':
                aln_type = atoi(optarg);
				if (aln_type < 0 || aln_type > 2)
                    return lsat_aln_usage();
                AP->aln_type = aln_type;
                break;

            case 'm': AP->match = atoi(optarg); break;
            case 'M': AP->mis = atoi(optarg); break;
            case 'O': AP->gapo = atoi(optarg); break;
            case 'E': AP->gape = atoi(optarg); break;

            case 'r': AP->res_mul_max = atoi(optarg); break;
            case 'V': AP->SV_len_thd = atoi(optarg); break;
            case 'g': AP->split_len = atoi(optarg);

            case 's': seed_program = atoi(optarg); break;
			case 'l': seed_len = atoi(optarg); break;
			case 'v': seed_inv = atoi(optarg); break;
			case 'p': per_aln = atoi(optarg); break;
			case 'd': min_thd = atoi(optarg); break;

            case 'o': strcpy(aln_result_f, optarg); break;

            case 'N': no_seed_aln = 1; break;
            case 'I': seed_info = 1; break;
            case 'A': strcpy(seed_result_f, optarg); break;	//seed alignment result break;

            default: return lsat_aln_usage();
        }
    }
    if (argc - optind != 2)
        return lsat_aln_usage();

    ref = strdup(argv[optind]);
    read =strdup(argv[optind+1]);
	if (seed_len != 0) {
		if (seed_inv == 0 || per_aln == 0 || min_thd == 0) {
			fprintf(stderr, "[lsat_aln] seed para error.\n");
			return lsat_aln_usage();
		}
		lsat_read_leve_n = 1;

		lsat_seed_len[0] = seed_len;
		lsat_seed_inv[0] = seed_inv;
		lsat_per_aln[0] = per_aln;
		lsat_min_thd[0] = min_thd;
	}

    lsat_aln_core(ref, read, seed_info, seed_program, no_seed_aln, seed_result_f, APP, AP);

    free(APP); free(AP);
    free(ref); free(read);
    return 0;
}
