#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>
#include <ctype.h>
#include "lsat_aln.h"
#include "bntseq.h"
#include "frag_check.h"
#include "split_mapping.h"
#include "kseq.h"
#include "kstring.h"
#include "./lsat_sam_parse/sam.h"

KSEQ_INIT(gzFile, gzread)
    
//                      0123456789
//#define FRAG_CON_STR "MPXLIDCRUS"
//#define F_MATCH 0
//#define F_SPLIT_MATCH 1
//#define F_MISMATCH 2
//#define F_LONG_MISMATCH 3
//#define F_INSERT 4
//#define F_DELETE 5
//#define F_CHR_DIF 6
//#define F_REVERSE 7
//#define F_UNCONNECT 8
//#define F_UNMATCH 9

int f_BCC_score_table[10] = {
      1,  //F_MATCH
      1,  //F_SPLIT_MATCH
      1,  //F_MISMATCH
     -1,  //F_LONG_MISMATCH
     -3,  //F_INSERT
     -3,  //F_DELETE
     -3,  //F_CHR_DIF
     -3,  //F_REVERSE
     -6,  //F_UNCONNECT
     -6,  //F_UNMATCH
};

int usage(void )		//aln usage
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   lsat aln [options] <ref_prefix> <in.fa/fq>\n\n");
    fprintf(stderr, "Options: -n            Do NOT excute soap2-dp program, when soap2-dp result existed already.\n");
    fprintf(stderr, "         -s            Seed information file has already existed.\n");
    fprintf(stderr, "         -a [STR]      The soap2-dp alignment result. When '-n' is used. [Def=\"seed_prefix.out.0\"]\n");
    fprintf(stderr, "         -m [INT][STR] The soap2-dp option. Maximun #errors allowed. [Def=3e]\n");
    fprintf(stderr, "         -l [INT]      The length of seed. [Def=100]\n");
    fprintf(stderr, "         -t [INT]      The alignment result type, [Def=1]\n");
    fprintf(stderr, "                           1 : Best coverage alignments.\n");
    fprintf(stderr, "                           2 : All valid alignments.\n");
    //fprintf(stderr, "         -o [STR]      The output file (SAM format). [Def=\"prefix_out.sam\"]\n");
    fprintf(stderr, "\n");
    return 1;
}

seed_msg *seed_init_msg(void)
{
    seed_msg *msg = (seed_msg*)malloc(sizeof(seed_msg));

    msg->n_read = 0;
    msg->read_m = READ_INIT_MAX;
    msg->n_seed = (int *)calloc(READ_INIT_MAX, sizeof(int));
    msg->read_name = (char **)malloc(READ_INIT_MAX * sizeof(char*));
    int i;
    for (i = 0; i < READ_INIT_MAX; ++i)
        msg->read_name[i] = (char*)malloc(1024 * sizeof(char));
    msg->last_len = (int *)malloc(READ_INIT_MAX * sizeof(int));
    msg->read_len = (int *)malloc(READ_INIT_MAX * sizeof(int));
    msg->seed_max = 0;
    msg->read_max_len = 0;

    return msg;
}

void seed_free_msg(seed_msg *msg)
{
    int i;
    for (i = 0; i < msg->read_m; ++i) free(msg->read_name[i]);
    free(msg->read_name);
    free(msg->n_seed);
    free(msg->last_len);
    free(msg->read_len);
    free(msg);
}

int split_seed(const char *prefix, seed_msg *s_msg, int seed_len)
{
    gzFile infp;
    kseq_t *seq;
    char out_f[1024], seed_head[1024], seed_seq[1024], seed_info[1024];
    FILE *outfp, *infofp;
    int m_read, n_seed, i;
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
    fprintf(infofp, "%d\n", seed_len);

    fprintf(stderr, "[lsat_aln] Spliting seed ... ");
    seed_seq[seed_len] = '\n';
    m_read = s_msg->read_m;
    while (kseq_read(seq) >= 0)
    {
        n_seed = ((seq->seq.l / seed_len) + 1) >> 1;    //XXX
        if (n_seed > s_msg->seed_max) s_msg->seed_max = n_seed;
        if (seq->seq.l > s_msg->read_max_len) s_msg->read_max_len = seq->seq.l; //XXX
        if (s_msg->n_read == m_read-1)
        {
            m_read <<= 1;
            if ((new_p = (int*)realloc(s_msg->n_seed, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->n_seed);
                fprintf(stderr, "[lsat_aln] Can't allocate more memory for n_seed[].\n");
                exit(1);
            }
            s_msg->n_seed = new_p;
            if ((new_p = (int*)realloc(s_msg->last_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->last_len);
                fprintf(stderr, "[lsat_aln] Can't allocate more memory for last_len[].\n");
                exit(-1);
            }
            s_msg->last_len = new_p;
            if ((new_p = (int*)realloc(s_msg->read_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_len);
                fprintf(stderr, "[lsat_aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_len = new_p;
            if ((new_p = (char**)realloc(s_msg->read_name, m_read * sizeof(char*))) == NULL)
            {
                free(s_msg->read_name);
                fprintf(stderr, "[lsat_aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_name = new_p;
            int i;
            for (i = m_read>>1; i < m_read; ++i)
                s_msg->read_name[i] = (char*)malloc(1024 * sizeof(char));
        }
        ++s_msg->n_read;
        s_msg->n_seed[s_msg->n_read] = s_msg->n_seed[s_msg->n_read-1] + n_seed;
        s_msg->last_len[s_msg->n_read] = seq->seq.l - (n_seed * 2 - 1) * seed_len;
        s_msg->read_len[s_msg->n_read] = seq->seq.l;
        strcpy(s_msg->read_name[s_msg->n_read], seq->name.s);
        s_msg->read_m = m_read;

        for (i = 0; i < n_seed; ++i)
        {
            sprintf(seed_head, ">%s_%d:%d\n", seq->name.s, i, i*seed_len*2);
            strncpy(seed_seq, seq->seq.s+i*seed_len*2, seed_len);
            seed_seq[seed_len+1] = '\0';
            fputs(seed_head, outfp);
            fputs(seed_seq, outfp);
        }
        fprintf(infofp, "%s %d %d %d\n", seq->name.s, n_seed, s_msg->last_len[s_msg->n_read], (int)seq->seq.l);
    }



    fprintf(stderr, "done.\n");
    gzclose(infp);
    fclose(outfp);
    fclose(infofp);
    kseq_destroy(seq);

    return 0;
}

int split_seed_info(const char *prefix, seed_msg *s_msg, int *seed_len)
{
    char seed_info[1024];
    char read_name[1024];
    FILE *infofp;
    int m_read, n_seed, last_len, len, n;
    void *new_p;

    strcpy(seed_info, prefix); strcat(seed_info, ".seed.info");
    if ((infofp = fopen(seed_info, "r")) == NULL)
    {
        fprintf(stderr, "[split seed] Can't open %s.\n", seed_info); 
        exit(-1);
    }
    m_read = s_msg->read_m;
    fprintf(stderr, "[last_aln] Parsing seeds' information ... ");
    if (fscanf(infofp, "%d", seed_len) == EOF)
    {
        fprintf(stderr, "[split seed] INFO file error.[1]\n");
        exit(-1);
    }
    while ((n = fscanf(infofp, "%s %d %d %d", read_name, &n_seed, &last_len, &len)) != EOF)
    {
        if (n != 4)
        {
            fprintf(stderr, "[split seed] INFO file error.[2]\n");
            exit(-1);
        }
        if (n_seed > s_msg->seed_max) s_msg->seed_max = n_seed;
        if (len > s_msg->read_max_len) s_msg->read_max_len = len;
        if (s_msg->n_read == m_read-1)
        {
            m_read <<= 1;
            if ((new_p = (int*)realloc(s_msg->n_seed, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->n_seed);
                fprintf(stderr, "[lsat aln] Can't allocate more memory for n_seed[].\n");
                exit(-1);
            }
            s_msg->n_seed =new_p;
            if ((new_p = (int*)realloc(s_msg->last_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->last_len);
                fprintf(stderr, "[lsat aln] Can't allocate more memory for last_len[].\n");
                exit(-1);
            }
            s_msg->last_len = new_p;
            if ((new_p = (int*)realloc(s_msg->read_len, m_read * sizeof(int))) == NULL)
            {
                free(s_msg->read_len);
                fprintf(stderr, "[lsat aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_len = new_p;
            //read_name
            if ((new_p = (char**)realloc(s_msg->read_name, m_read * sizeof(char*))) == NULL)
            {
                free(s_msg->read_name);
                fprintf(stderr, "[lsat_aln] Can't allocate more memory for read_len[].\n");
                exit(-1);
            }
            s_msg->read_name = new_p;
            int i;
            for (i = m_read>>1; i < m_read; ++i)
                s_msg->read_name[i] = (char*)malloc(1024 * sizeof(char));
        }
        ++s_msg->n_read;
        s_msg->n_seed[s_msg->n_read] = s_msg->n_seed[s_msg->n_read-1] + n_seed;
        s_msg->last_len[s_msg->n_read] = last_len;
        s_msg->read_len[s_msg->n_read] = len;
        strcpy(s_msg->read_name[s_msg->n_read], read_name);
        s_msg->read_m = m_read;
        if (last_len != len - (n_seed * 2 - 1) * (*seed_len))
        {
            fprintf(stderr, "[split seed] INFO file error.[3]\n");
            exit(-1);
        }
    }

    fprintf(stderr, "done.\n");
    fclose(infofp);
    return 0;
}

sam_msg *sam_init_msg(void)
{
    sam_msg *s = (sam_msg*)malloc(sizeof(sam_msg));
    s->sam_n = 0;
    s->sam_m = 1;
    s->sam = (sam_t*)malloc(sizeof(sam_t));
    s->sam->cigar_s = (kstring_t*)calloc(1, sizeof(kstring_t));
    return s;
}

aln_msg *aln_init_msg(int seed_max)
{
    aln_msg *msg;
    int i,j;
    msg = (aln_msg*)malloc(seed_max * sizeof(aln_msg));
    for (i = 0; i < seed_max; ++i)		//drop away seed whose number of alignments > PER_ALN_N
    {
        msg[i].read_id = -1;
        msg[i].n_aln = 0;
        msg[i].skip = 0;
        msg[i].at = (aln_t*)malloc(PER_ALN_N * sizeof(aln_t));
        for (j = 0; j < PER_ALN_N; ++j)
        {
            msg[i].at[j].cigar = (uint32_t*)malloc(7 * sizeof(uint32_t));//XXX default value for 3-ed
            msg[i].at[j].cigar_len = 0;
            msg[i].at[j].cmax = 7;
            msg[i].at[j].bmax = 0;
        }
    }
    return msg;
}

void aln_free_msg(aln_msg *a_msg, int seed_max)	//a_msg[seed_max]
{
    int i,j;
    for (i = 0; i < seed_max; ++i)
    {
        for (j = 0; j < PER_ALN_N; ++j)
        {
            free(a_msg[i].at[j].cigar);
        }
        free(a_msg[i].at);
    }
    free(a_msg);
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
    for (s = s_cigar; *s; )
    {
        x = strtol(s, &t, 10);	
        if (x == 0)
        {
            fprintf(stderr, "%s\n",s);
            fprintf(stderr, "[lsat_aln] Cigar ERROR 1.\n");
            exit(-1);
        }
        op = toupper(*t);
        switch (op)
        {
            case 'M':	op = CMATCH;	break;
            case 'I':	op = CINS;		bi += x;	break;
            case 'D':	op = CDEL;		bd += x;	break;
            case 'N':	op = CREF_SKIP;		bd += x;	break;
            case 'S':	op = CSOFT_CLIP;		bi += x;	break;
            case 'H':	op = CHARD_CLIP;		bd += x;	break;
            case 'P':	op = CPAD;		bd += x;	break;
            case '=':	op = CEQUAL;	break;
            case 'X':	op = CDIFF;	break;
            case 'B':	op = CBACK;		bi += x;	break;	
            default:	fprintf(stderr, "[lsat_aln] Cigar ERROR 2.\n"); exit(-1); break;
        }
        if (a_msg[seed_i].at[aln_i].cigar_len == a_msg[seed_i].at[aln_i].cmax)
        {
            a_msg[seed_i].at[aln_i].cmax <<= 2 ;
            a_msg[seed_i].at[aln_i].cigar = (uint32_t*)realloc(a_msg[seed_i].at[aln_i].cigar, a_msg[seed_i].at[aln_i].cmax * sizeof(uint32_t));
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
        int read_id, int chr, int64_t offset, 
        char srand, /*int edit_dis,*/ char *cigar)
{   //read_x: (除去unmap和repeat的)read序号, aln_y: read对应的比对结果序号(从1开始)
    if (aln_y > PER_ALN_N) {
        fprintf(stderr, "[lsat_aln] setAmsg ERROR!\n");
        exit(0);
    }
    //printf("%d %d: %d %d %lld %c %s\n", read_x, aln_y, read_id, chr, offset, srand, cigar);
    a_msg[read_x-1].read_id = read_id;			//from 1
    a_msg[read_x-1].at[aln_y-1].chr = chr;
    a_msg[read_x-1].at[aln_y-1].offset = offset;	//1-base
    a_msg[read_x-1].at[aln_y-1].nsrand = ((srand=='+')?1:-1);
    //a_msg[read_x-1].at[aln_y-1].edit_dis = edit_dis;
    a_msg[read_x-1].n_aln = aln_y;
    setCigar(a_msg, read_x-1, aln_y-1,  cigar);
}

//for best coverage and connect
//add "overlap ins"
int get_abs_dis(aln_msg *a_msg, int pre, int pre_a, int i, int j, int *flag, int seed_len)    //(i,j)对应节点，来自pre的第pre_a个aln
{
    if (pre == -1 || i == -1) { *flag = F_MATCH; return 0; }	//for bound node
    if (pre == i)// XXX overlap???
    {
        if (pre_a == j)
            *flag = F_MATCH; 
        else
            *flag = F_UNCONNECT;
        return 0;
    }
    if (a_msg[i].at[j].chr != a_msg[pre].at[pre_a].chr || a_msg[i].at[j].nsrand != a_msg[pre].at[pre_a].nsrand)	//different chr or different srnad
    {
        *flag = F_CHR_DIF;
        return 0;//PRICE_DIF_CHR;
    }
    int64_t exp = a_msg[pre].at[pre_a].offset + a_msg[pre].at[pre_a].nsrand * (a_msg[i].read_id - a_msg[pre].read_id) * 2 * seed_len;	
    int64_t act = a_msg[i].at[j].offset;
    int64_t dis = a_msg[pre].at[pre_a].nsrand * ((a_msg[pre].read_id < a_msg[i].read_id)?(act-exp):(exp-act)) - (((a_msg[pre].at[pre_a].nsrand) * (a_msg[pre].read_id-a_msg[i].read_id) < 0)?(a_msg[pre].at[pre_a].len_dif):(a_msg[i].at[j].len_dif));

    if (dis <= 10 && dis >= -10) {
        if (abs(a_msg[pre].read_id - a_msg[i].read_id) == 1) *flag = F_MATCH;
        else if (abs(a_msg[pre].read_id - a_msg[i].read_id) < 10) *flag = F_MISMATCH; // XXX long dis
        else *flag = F_LONG_MISMATCH;
    } else if (dis > 10 && dis < DEL_THD) *flag = F_DELETE;
    else if (dis < -10 && dis >= (0-(abs(a_msg[i].read_id-a_msg[pre].read_id)*2-1)*seed_len)) *flag = F_INSERT;
    //else if (dis < -10 && dis >= (0-(abs(a_msg[i].read_id-a_msg[pre].read_id)*2+1)*seed_len)) *flag = F_INSERT; // overlap ins
    else { *flag = F_UNCONNECT; return 0; }
    //printf("getdis: readid: %d %d dis %d flag %d\n", a_msg[pre].read_id, a_msg[i].read_id, dis, *flag);
    dis=abs(dis);
    dis += (((a_msg[pre].at[pre_a].nsrand) * (a_msg[pre].read_id - a_msg[i].read_id) < 0)? a_msg[i].at[j].cigar_len : a_msg[pre].at[pre_a].cigar_len);
    return dis; // XXX dis..
}

void new_line(int seed_i, int **line, int *path_end, int path_n)
{
    line[path_n][0] = seed_i;
    path_end[path_n] = 1;
}

void frag_new_line(int seed_i, int aln_i, line_node **line, int *line_end, int path_n)
{
    line[path_n][0].x = seed_i;
    line[path_n][0].y = aln_i;
    line_end[path_n] = 1;
}

void copy_line(int **line, int from, int to, int *path_end)
{
    int i;
    for (i = 0; i < path_end[from]; ++i)
        line[to][path_end[to]+i] = line[from][i];
    path_end[to] += path_end[from];
}

int frag_add_seed(aln_msg *a_msg, int seed_i, int aln_i, 
        int *path_n, line_node **line, int *line_end, 
        int seed_len)
{
    int con_flag;
    if (*path_n > 0)
    {
        if (seed_i == line[*path_n-1][line_end[*path_n-1]-1].x + 1)
        {
            get_abs_dis(a_msg, line[*path_n-1][line_end[*path_n-1]-1].x, line[*path_n-1][line_end[*path_n-1]-1].y, seed_i, aln_i, &con_flag, seed_len);	
            if (con_flag == F_MATCH)
            {
                line[*path_n-1][line_end[*path_n-1]].x = seed_i;
                line[*path_n-1][line_end[*path_n-1]].y = aln_i;
                ++line_end[*path_n-1];
                return 0;
            }
        }
    }
    line[*path_n][0].x = seed_i;
    line[*path_n][0].y = aln_i;
    line_end[*path_n] = 1;
    ++(*path_n);
    return 0;
}

int frag_DP_init(frag_DP_node *f_node, 
        line_node **line, int *line_end, 
        int path_i, line_node from, 
        aln_msg *a_msg, int seed_len)
{
    f_node[path_i].seed_num = line_end[path_i];
    f_node[path_i].n_line = 1;
    f_node[path_i].line_i[0] = path_i;
    f_node[path_i].from_i = -1;
    int con_flag;
    get_abs_dis(a_msg, from.x, from.y, line[path_i][0].x, line[path_i][0].y, &con_flag, seed_len);
    f_node[path_i].score = (con_flag == F_MATCH ? (line_end[path_i]):(line_end[path_i] - F_SV_PEN));
    return 0;
}

int frag_update_line(frag_DP_node *f_node, 
        line_node **line, int *line_end, 
        int path_i, aln_msg *a_msg, int seed_len)
{
    int i, pre, pre_i, cur, cur_i, con_flag;
    int max_score = f_node[path_i].score, max_from = -1;
    for (i = 0; i < path_i; ++i)
    {
        if (line_end[i] == 0 ) continue;
        //1. check if line[i] is allowed to be combined with line[path_i]
        pre = line[i][line_end[i]-1].x;
        pre_i = line[i][line_end[i]-1].y;
        cur = line[path_i][0].x;
        cur_i = line[path_i][0].y;

        get_abs_dis(a_msg, pre, pre_i, cur, cur_i, &con_flag, seed_len);
        if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
            continue;
        //2. check if line[i] is the most suitable one
        if (con_flag == F_MATCH)
        {
            if (f_node[i].score + f_node[path_i].seed_num + F_MATCH_SCORE > max_score)
            {
                max_from = i;
                max_score = f_node[i].score + f_node[path_i].seed_num;
            }
        }
        else if (f_node[i].score + f_node[path_i].seed_num - F_SV_PEN > max_score)
        {
            max_score = f_node[i].score + f_node[path_i].seed_num - F_SV_PEN;
            max_from = i;
        }
    }
    if (max_from != -1)
    {
        f_node[path_i].from_i = max_from;
        f_node[path_i].score = max_score;
        f_node[path_i].seed_num += f_node[max_from].seed_num;
        for (i = 0; i < f_node[max_from].n_line; ++i)
            f_node[path_i].line_i[i] = f_node[max_from].line_i[i];
        f_node[path_i].line_i[i] = path_i;
        f_node[path_i].n_line = i+1;
    }

    return 0;
}

int frag_copy_line(line_node **line, int *line_end, int from, int to)
{
    int i;
    for (i = 0; i < line_end[from]; ++i)
    {
        line[to][line_end[to]+i] = line[from][i];
    }
    line_end[to] += line_end[from];
    return 0;
}

int mini_frag_main_line(aln_msg *a_msg, 
        line_node left, line_node right, 
        int seed_len, int n_seed, 
        line_node **line, int *line_end, 
        frag_DP_node *f_node)
{
    //two bounds
    int start = left.x+1, end = right.x-1;

    int i, j, con_flag;
    int path_n = 0;
    line_end[0] = 0;
    for (i = start; i <= end; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            if (left.x != -1)
            {
                get_abs_dis(a_msg, i, j, left.x, left.y, &con_flag, seed_len);
                if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
                    continue;
            }
            if (right.x != n_seed)
            {
                get_abs_dis(a_msg, i, j, right.x, right.y, &con_flag, seed_len);
                if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
                    continue;
            }
            //F_MATCH or SV
            frag_add_seed(a_msg, i, j, &path_n, line, line_end, seed_len);
        }
    }
    if (path_n == 0) return 0;
    //XXX bound
    line_end[path_n] = 1;
    line[path_n][0] = right;
    for (i = 0; i <= path_n; ++i)
        frag_DP_init(f_node, line, line_end, i, left, a_msg, seed_len);
    //XXX last node: right
    for (i = 1; i <= path_n; ++i)
        frag_update_line(f_node, line, line_end, i, a_msg, seed_len);
    line_end[path_n] = 0;
    for (i = 0; i < f_node[path_n].n_line-1; ++i)
        frag_copy_line(line, line_end, f_node[path_n].line_i[i], path_n);
    //printf("mini path:\t");
    //for (i = 0; i < line_end[path_n]; ++i)
    //	printf("%d %d\t", line[path_n][i].x, line[path_n][i].y);
    return path_n;
}

int frag_copy_main_line(line_node **line, int *line_end, int to, line_node **_line, int *_line_end, int from)
{
    int i, j;
    for (i = 0; i < line_end[to]; ++i)
    {
        if (line[to][i].x > _line[from][0].x)
            break;
    }
    //copy the tail of 'to' to 'to'
    for (j = line_end[to]-1+_line_end[from]; j > i + _line_end[from]-1; --j)
        line[to][j] = line[to][j-_line_end[from]];
    //copy 'from' to 'to'
    for (j = i; j < i + _line_end[from]; ++j)
        line[to][j] = _line[from][j-i];
    line_end[to] += _line_end[from];
    return 0;
}

int frag_main_line(aln_msg *a_msg, int n_seed, int seed_len, line_node **line, int *line_end, frag_DP_node *f_node)
{
    int path_n, i;
    int min_len;

    path_n = 0;
    line_end[0] = 0;
    min_len = a_msg[0].n_aln;
    //match-line of uniquely-aln seeds
    for (i = 0; i < n_seed; ++i)
    {
        if (a_msg[i].n_aln == 1)
            frag_add_seed(a_msg, i, 0, &path_n, line, line_end, seed_len);
        else if (min_len != 1 && a_msg[i].n_aln < min_len && a_msg[i].n_aln >= 2)
            min_len = a_msg[i].n_aln;
    }
    /*
    //match-line of min_len-aln seeds
    if (path_n == 0)
    {
    for (i = 0; i < n_seed; ++i)
    {
    if (a_msg[i].n_aln == min_len)
    {
    for (j = 0; j < min_len; ++j)
    {
    frag_add_seed(a_msg, i, j, &path_n, line, line_end, seed_len);	
    }
    }
    }
    }
    */
    /*//add multi-aln seeds to the unique(min_len)-aln match-line
      for (i = 0; i < n_seed; ++i)
      {
      if (a_msg[i].n_aln > min_len)
      {
    //j or k XXX
    for (j = 0; j < a_msg[i].n_aln; ++j)
    {
    for (k = path_n-1; k >=0; --k)
    {
    get_abs_dis(a_msg, line[i][line_end[i]-1].x,)
    }
    }
    }
    }*/

    for (i = 0; i < path_n; ++i)
        if (line_end[i] != 0) frag_DP_init(f_node, line, line_end, i, (line_node){-1,-1}, a_msg, seed_len);
    for (i = 1; i < path_n; ++i)
        if (line_end[i] != 0) frag_update_line(f_node, line, line_end, i, a_msg, seed_len);


    //find backtrack node
    int max_node = -1, max_score = 0;
    for (i = 0; i < path_n; ++i)
    {
        if (line_end[i] == 0) continue;
        if (f_node[i].score >= max_score)
        {
            max_node = i;
            max_score = f_node[i].score;
        }
    }
    line_end[path_n] = 0;
    if (max_node != -1)
    {
        for (i = 0; i < f_node[max_node].n_line; ++i)
            frag_copy_line(line, line_end, f_node[max_node].line_i[i], path_n);
    }

#ifdef __DEBUG__
    printf("line:\n");
    for (i = 0; i < line_end[path_n]; ++i)
    {
        printf("%d %d\n", line[path_n][i].x, line[path_n][i].y);
    }
#endif
    //for exist blank gap, use mini_main_line to fill it
    line_node **_line; int *_line_end; frag_DP_node *_f_node;
    int line_1, line_2, mini_i;

    //for _line: n_seed is enough?
    _line = (line_node**)malloc(n_seed * sizeof(line_node*));
    _line_end = (int*)malloc(n_seed * sizeof(int));
    _f_node = (frag_DP_node*)malloc(n_seed * sizeof(frag_DP_node));
    for (i = 0; i < n_seed; ++i)
    {
        _line[i] = (line_node*)malloc(n_seed * sizeof(line_node));
        _f_node[i].line_i = (int*)malloc(n_seed * sizeof(int));
    }

    //blank between left bound and first node-line
    line_2 = f_node[max_node].line_i[0];
    if (line[line_2][0].x > 0)	//blank exists
    {
        mini_i = mini_frag_main_line(a_msg, (line_node){-1,-1}, line[line_2][0], seed_len, n_seed, _line, _line_end, _f_node);
        if (_line_end[mini_i] != 0)
            frag_copy_main_line(line, line_end, path_n, _line, _line_end, mini_i);
    }
    //blank between line-nodes
    for (i = 0; i < f_node[max_node].n_line-1; ++i)
    {
        line_1 = f_node[max_node].line_i[i];
        line_2 = f_node[max_node].line_i[i+1];
        if (line[line_1][line_end[line_1]-1].x < line[line_2][line_end[line_2]-1].x-1)	//blank exists
        {
            mini_i = mini_frag_main_line(a_msg, line[line_1][line_end[line_1]-1], line[line_2][0], seed_len, n_seed, _line, _line_end, _f_node);
            if (_line_end[mini_i] != 0)
                frag_copy_main_line(line, line_end, path_n, _line, _line_end, mini_i);
        }
    }
    //blank betweed last node-line and right bound
    line_1 = f_node[max_node].line_i[i];
    if (line[line_1][line_end[line_1]-1].x < n_seed-1)	//blank exists
    {
        mini_i = mini_frag_main_line(a_msg, line[line_1][line_end[line_1]-1], (line_node){n_seed,0}, seed_len, n_seed, _line, _line_end, _f_node);
        if (_line_end[mini_i] != 0)
            frag_copy_main_line(line, line_end, path_n, _line, _line_end, mini_i);
    }

    //free variable
    for (i = 0; i < n_seed; ++i) { free(_line[i]); free(_f_node[i].line_i); }
    free(_line); free(_line_end); free(_f_node);

    return path_n;
}

void f_node_set(frag_dp_node ***f_node, int f_x, int f_y, 
        line_node from, int seed_i, int aln_i, 
        int score, int dis_pen,
        uint8_t match_flag, int dp_flag, int backtrack_flag, 
        int node_n, int peak_value, line_node peak, 
        int next_trigger_n, int pre_trigger_n)
{
    (*f_node)[f_x][f_y].from = from;
    (*f_node)[f_x][f_y].seed_i = seed_i;
    (*f_node)[f_x][f_y].aln_i = aln_i;

    (*f_node)[f_x][f_y].score = score;
    (*f_node)[f_x][f_y].dis_pen = dis_pen;

    (*f_node)[f_x][f_y].match_flag = match_flag;
    (*f_node)[f_x][f_y].dp_flag = dp_flag;
    (*f_node)[f_x][f_y].backtrack_flag = backtrack_flag;

    (*f_node)[f_x][f_y].node_n = node_n;

    (*f_node)[f_x][f_y].peak_value = peak_value;
    (*f_node)[f_x][f_y].peak = peak;

    (*f_node)[f_x][f_y].next_trigger_n = next_trigger_n;
    (*f_node)[f_x][f_y].pre_trigger_n = pre_trigger_n;
}

int frag_dp_init(frag_dp_node **f_node, 
                 aln_msg *a_msg, int seed_i, 
                 line_node from, int seed_len, int dp_flag)
{
    int i;
    if (from.x == -1)	//UNLIMITED
    {
        for (i = 0; i < a_msg[seed_i].n_aln; ++i)
        {
            //init-score: 0 == 1 + 1 - F_SV_PEN
            //f_node[seed_i][i] = (frag_dp_node){from, seed_i, i, 0, F_MATCH, dp_flag, 1, 1, 0, (line_node){-1,-1}, 0, 0};
            f_node_set(&f_node, seed_i, i, from, seed_i, i, 0, a_msg[seed_i].at[i].cigar_len, F_MATCH, dp_flag, 1, 1, 0, (line_node){-1, -1}, 0, 0);
        }
    }
    else
    {
        int con_flag, con_score, dis;
        for (i = 0; i < a_msg[seed_i].n_aln; ++i)
        {
            
            dis = get_abs_dis(a_msg, from.x, from.y, seed_i, i, &con_flag, seed_len);
            //											XXX MATCH, MISMATCH: same score?
            //											mismatch too long, socre??? XXX
            con_score = f_BCC_score_table[con_flag];

            //f_node_set(&f_node, seed_i, i, from, seed_i, i, 1+con_score, dis, con_flag, dp_flag, 1, 1, 0, (line_node){-1,-1}, 0, 0);

            if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
                f_node_set(&f_node, seed_i, i, from, seed_i, i, 1+con_score, dis, con_flag, dp_flag, 1, 1, 0, (line_node){-1,-1}, 0, 0);
            else
                f_node_set(&f_node, seed_i, i, (line_node){-1,0}, seed_i, i, -1, 0, con_flag, 0-dp_flag, 1, 0, 0, (line_node){-1,-1}, 0, 0);
        }
    }
    return 0;
}

int frag_dp_multi_init(frag_dp_node **f_node, 
        aln_msg *a_msg, 
        int seed_i, line_node from, 
        int seed_len, int dp_flag)
{
    int i;
    if (from.x == -1)	//UNLIMITED
    {
        for (i = 0; i < a_msg[seed_i].n_aln; ++i)
        {
            //init-score: 0 == 1 + 1 - F_SV_PEN
            //f_node[seed_i][i] = (frag_dp_node){from, seed_i, i, 0, F_MATCH, dp_flag, 1, 1, 0, (line_node){-1,-1}, 0, 0};
            f_node_set(&f_node, seed_i, i, from, seed_i, i, 0, 0, F_MATCH, dp_flag, 1, 1, 0, (line_node){-1,-1}, 0, 0);
        }
    }
    else
    {
        int con_flag, dis;
        for (i = 0; i < a_msg[seed_i].n_aln; ++i)
        {
            dis = get_abs_dis(a_msg, from.x , from.y, seed_i, i, &con_flag, seed_len);
            if (con_flag == F_MATCH || con_flag == F_MISMATCH)
                f_node_set(&f_node, seed_i, i, from, seed_i, i, 2+F_MATCH_SCORE, dis, con_flag, dp_flag, 1, 1, 0, (line_node){-1,-1}, 0, 0);
            else
                f_node_set(&f_node, seed_i, i, (line_node){-1,0}, seed_i, i, -1, 0, con_flag, 0-dp_flag, 1, 0, 0, (line_node){-1,-1}, 0, 0);
        }
    }
    return 0;
}

//pruning	XXX
/*int frag_dp_update(frag_dp_node **f_node, 
        aln_msg *a_msg, 
        int seed_i, int aln_i, 
        int start, int seed_len, int dp_flag)
{
    int i, j, con_flag, dis;
    line_node max_from;
    int max_score, max_flag, max_dis;

    max_from = f_node[seed_i][aln_i].from;
    max_score = f_node[seed_i][aln_i].score;
    max_dis = 0;

    for (i = seed_i - 1; i >= start; --i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            if (f_node[i][j].dp_flag == dp_flag)
            {
                dis = get_abs_dis(a_msg, i, j, seed_i, aln_i, &con_flag, seed_len);
                if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
                    continue;
                //											XXX MATCH, MISMATCH: same score?
                if ((f_node[i][j].score + 1 + ((con_flag == F_MATCH||con_flag == F_MISMATCH)?F_MATCH_SCORE:0-F_SV_PEN)) > max_score)	
                {
                    max_from = (line_node){i,j};
                    max_score = f_node[i][j].score + 1 + (con_flag==F_MATCH||con_flag ==F_MISMATCH?F_MATCH_SCORE:0-F_SV_PEN);
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                }
                else if ((f_node[i][j].score + 1 + ((con_flag == F_MATCH||con_flag == F_MISMATCH)?F_MATCH_SCORE:0-F_SV_PEN)) == max_score)
                {
                    if (f_node[i][j].dis_pen + dis < max_dis)
                    {
                        max_from = (line_node){i,j};
                        max_score = f_node[i][j].score + 1 + (con_flag==F_MATCH||con_flag ==F_MISMATCH?F_MATCH_SCORE:0-F_SV_PEN);
                        max_flag = con_flag;
                        max_dis = f_node[i][j].dis_pen + dis;
                    }
                }
            }
        }
    }
    if (max_from.x != f_node[seed_i][aln_i].from.x || max_from.y != f_node[seed_i][aln_i].from.y)
    {
        f_node[seed_i][aln_i].from = max_from;
        f_node[seed_i][aln_i].score = max_score;
        f_node[seed_i][aln_i].dis_pen = max_dis;
        f_node[seed_i][aln_i].match_flag = max_flag;
        f_node[seed_i][aln_i].node_n += f_node[max_from.x][max_from.y].node_n;

        if (f_node[max_from.x][max_from.y].backtrack_flag == 1) //local maximum value
        {
            if (max_score >= f_node[max_from.x][max_from.y].score)
            {
                f_node[max_from.x][max_from.y].backtrack_flag = 0;
            }
            else
            {
                f_node[seed_i][aln_i].backtrack_flag = 0;
                f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].score;
                f_node[seed_i][aln_i].peak = max_from;
            }
        }
        else    
        {
            if (f_node[max_from.x][max_from.y].peak.x != -1) //there exist a peak node
            {
                if (max_score >= f_node[max_from.x][max_from.y].peak_value)
                {
                    f_node[f_node[max_from.x][max_from.y].peak.x][f_node[max_from.x][max_from.y].peak.y].backtrack_flag = 0;
                }
                else
                {
                    f_node[seed_i][aln_i].backtrack_flag = 0;
                    f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].peak_value;
                    f_node[seed_i][aln_i].peak = f_node[max_from.x][max_from.y].peak;
                }
            }
            else // NO peak node
            {
                if (max_score < f_node[max_from.x][max_from.y].score)
                {
                    f_node[seed_i][aln_i].backtrack_flag = 0;
                    f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].score;
                    f_node[seed_i][aln_i].peak = max_from;
                }
            }
        }
    }
    return 0;
}*/
//for best coverage and connect
//    similar to MEM, match/mismatch have the highest priority in dp.
int frag_dp_update(frag_dp_node **f_node, 
                   aln_msg *a_msg, 
                   int seed_i, int aln_i, 
                   int start, int seed_len, int dp_flag)
{
    int i, j, con_flag, con_score, dis;
    line_node max_from;
    int max_score, max_flag, max_dis;

    max_from = f_node[seed_i][aln_i].from;
    max_score = f_node[seed_i][aln_i].score;
    max_dis = 0;
    max_flag = 0;

    for (i = seed_i - 1; i >= start; --i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            if (f_node[i][j].dp_flag == dp_flag)
            {
                dis = get_abs_dis(a_msg, i, j, seed_i, aln_i, &con_flag, seed_len);
                if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
                    continue;
                //											XXX MATCH, MISMATCH: same score?
                con_score = f_BCC_score_table[con_flag];

                if (con_flag == F_MATCH || con_flag == F_MISMATCH) // MEM // XXX cigar-diff, when score equal with each other
                {
                    max_from = (line_node){i,j};
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                    goto UPDATE;
                }

                if (f_node[i][j].score + 1 + con_score > max_score)	
                {
                    max_from = (line_node){i,j};
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                }
                else if (f_node[i][j].score + 1 + con_score == max_score)
                {
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
    if (max_from.x != f_node[seed_i][aln_i].from.x || max_from.y != f_node[seed_i][aln_i].from.y)
    {
        f_node[seed_i][aln_i].from = max_from;
        f_node[seed_i][aln_i].score = max_score;
        f_node[seed_i][aln_i].dis_pen = max_dis;
        f_node[seed_i][aln_i].match_flag = max_flag;
        f_node[seed_i][aln_i].node_n += f_node[max_from.x][max_from.y].node_n;

        if (f_node[max_from.x][max_from.y].backtrack_flag == 1) //local maximum value
        {
            if (max_score >= f_node[max_from.x][max_from.y].score)
            {
                f_node[max_from.x][max_from.y].backtrack_flag = 0;
            }
            else
            {
                f_node[seed_i][aln_i].backtrack_flag = 0;
                f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].score;
                f_node[seed_i][aln_i].peak = max_from;
            }
        }
        else    
        {
            if (f_node[max_from.x][max_from.y].peak.x != -1) //there exist a peak node
            {
                if (max_score >= f_node[max_from.x][max_from.y].peak_value)
                {
                    f_node[f_node[max_from.x][max_from.y].peak.x][f_node[max_from.x][max_from.y].peak.y].backtrack_flag = 0;
                }
                else
                {
                    f_node[seed_i][aln_i].backtrack_flag = 0;
                    f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].peak_value;
                    f_node[seed_i][aln_i].peak = f_node[max_from.x][max_from.y].peak;
                }
            }
            else // NO peak node
            {
                if (max_score < f_node[max_from.x][max_from.y].score)
                {
                    f_node[seed_i][aln_i].backtrack_flag = 0;
                    f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].score;
                    f_node[seed_i][aln_i].peak = max_from;
                }
            }
        }
    }
    return 0;
}

int frag_dp_multi_update(frag_dp_node **f_node, 
        aln_msg *a_msg, 
        int seed_i, int aln_i, 
        int start, int seed_len, int dp_flag)
{
    int i, j, con_flag, dis;
    line_node max_from;
    int max_score, max_flag, max_dis; // TODO: store as a struct or sth.

    max_from = f_node[seed_i][aln_i].from;
    max_score = f_node[seed_i][aln_i].score;
    max_dis = 0;

    for (i = seed_i - 1; i >= start; --i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            if (f_node[i][j].dp_flag == dp_flag)
            {
                dis = get_abs_dis(a_msg, i, j, seed_i, aln_i, &con_flag, seed_len);
                if (con_flag != F_MISMATCH && con_flag != F_MATCH)
                    continue;
                //											XXX MATCH, MISMATCH: same score?
                if ((f_node[i][j].score + 1 + F_MATCH_SCORE) > max_score)	
                {
                    max_from = (line_node){i,j};
                    max_score = f_node[i][j].score + 1 + F_MATCH_SCORE;
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                }
                else if ((f_node[i][j].score + 1 + F_MATCH_SCORE) == max_score)
                {
                    if (f_node[i][j].dis_pen + dis < max_dis)
                    {
                        max_from = (line_node){i,j};
                        max_score = f_node[i][j].score + 1 + F_MATCH_SCORE;
                        max_flag = con_flag;
                        max_dis = f_node[i][j].dis_pen + dis;
                    }
                }
            }
        }
    }
    if (max_from.x != f_node[seed_i][aln_i].from.x || max_from.y != f_node[seed_i][aln_i].from.y)
    {
        f_node[seed_i][aln_i].from = max_from;
        f_node[seed_i][aln_i].score = max_score;
        f_node[seed_i][aln_i].dis_pen = max_dis;
        f_node[seed_i][aln_i].match_flag = max_flag;
        f_node[seed_i][aln_i].node_n += f_node[max_from.x][max_from.y].node_n;

        if (f_node[max_from.x][max_from.y].backtrack_flag == 1) //local maximum value
        {
            if (max_score >= f_node[max_from.x][max_from.y].score)
            {
                f_node[max_from.x][max_from.y].backtrack_flag = 0;
            }
            else
            {
                f_node[seed_i][aln_i].backtrack_flag = 0;
                f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].score;
                f_node[seed_i][aln_i].peak = max_from;
            }
        }
        else    
        {
            if (f_node[max_from.x][max_from.y].peak.x != -1) //there exist a peak node
            {
                if (max_score >= f_node[max_from.x][max_from.y].peak_value)
                {
                    f_node[f_node[max_from.x][max_from.y].peak.x][f_node[max_from.x][max_from.y].peak.y].backtrack_flag = 0;
                }
                else
                {
                    f_node[seed_i][aln_i].backtrack_flag = 0;
                    f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].peak_value;
                    f_node[seed_i][aln_i].peak = f_node[max_from.x][max_from.y].peak;
                }
            }
            else // NO peak node
            {
                if (max_score < f_node[max_from.x][max_from.y].score)
                {
                    f_node[seed_i][aln_i].backtrack_flag = 0;
                    f_node[seed_i][aln_i].peak_value = f_node[max_from.x][max_from.y].score;
                    f_node[seed_i][aln_i].peak = max_from;
                }
            }
        }
    }
    return 0;
}
//for multi-dp-lines
//line, line_end: start at 1
//left and right are 'MIN_FLAG'	XXX
int frag_mini_dp_multi_line(frag_dp_node **f_node, 
        aln_msg *a_msg, int seed_len, 
        line_node left, line_node right, 
        line_node **line, int *line_end, int max_multi,
        int _head,  int _tail, int n_seed)
{
    //line_node head = ((_head)?left:(line_node){-1,0});
    line_node head = (line_node){-1, 0};
    int start, end;
    int i, j, k, l, mini_dp_flag = MULTI_FLAG;
    int l_i, max_score, node_i;
    line_node _right, _left;
    int candi_n; line_node candi[max_multi];

    if (_head) start = left.x; else start = left.x+1;
    if (_tail) end = right.x; else end = right.x-1;

    //first dp int
    for (i = start; i <= end; ++i)
        frag_dp_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);
    //dp update
    for (i = start+1; i <= end; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            frag_dp_update(f_node, a_msg, i, j, start, seed_len, mini_dp_flag);
        }
    }

    //store and sort candidate node, the biggest 'max-multi' cand-nodes
    candi_n=0;
    for (i = right.x-1; i > left.x; --i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            //candi-nodes
            if (f_node[i][j].backtrack_flag == 1 && f_node[i][j].dp_flag == mini_dp_flag && f_node[i][j].score > 0)
            {
                //sort candi-nodes in 'candi', big->small
                for (k = 0; k < candi_n; ++k)
                {
                    if (f_node[i][j].score > f_node[candi[k].x][candi[k].y].score)
                        break;
                }
                if (k >= max_multi) continue;
                if (candi_n < max_multi) candi_n++;
                for (l = candi_n-1; l > k; --l)
                {
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
    if (_tail)
    {
        //right.x, right.y -> left.x, left.y
        if (_head)
        {
            node_i = f_node[right.x][right.y].node_n - 3;
            if (node_i >= 0) 
            {
                line_end[0] = node_i+1;
                _right = right;
                _left = f_node[_right.x][_right.x].from;
                while (_left.x != head.x && _left.x != left.x)
                {
                    if (node_i < 0) { fprintf(stderr, "[frag mini multi dp] node_i < 0, BUG.\n"); exit(0); }
                    line[0][node_i] = _left;
                    --node_i;
                    _right = _left;
                    _left = f_node[_right.x][_right.y].from;
                }
                if (node_i >= 0) { fprintf(stderr, "[frag mini multi dp] node_i >= 0, BUG.\n"); exit(0); }
                if (_right.x != left.x ||_right.y != left.y)		//'anchor-match' line
                    line_end[0] = 0;
            }
        }
        //right.x, right.y
        else
        {
            node_i = f_node[right.x][right.y].node_n - 2;
            if (node_i >= 0)
            {
                line_end[0] = node_i+1;
                _right = right;
                _left = f_node[_right.x][_right.x].from;
                while (_left.x != head.x)
                {
                    if (node_i < 0) { fprintf(stderr, "[frag mini multi dp] node_i < 0, BUG.\n"); exit(0); }
                    line[0][node_i] = _left;
                    --node_i;
                    _right = _left;
                    _left = f_node[_right.x][_right.y].from;
                }
                if (node_i >= 0) { fprintf(stderr, "[frag mini multi dp] node_i >= 0, BUG.\n"); exit(0); }
            }
        }
    }
    else if (_head)
    {
        i = 0;
        while (i < candi_n)
        {
            node_i = f_node[candi[i].x][candi[i].y].node_n - 2;
            if (node_i >= 0)
            {
                line_end[0] = node_i+1;
                _right = candi[i];
                _left = f_node[_right.x][_right.y].from;
                while (_left.x != head.x && _left.x != left.x)
                {
                    if (node_i < 0) { fprintf(stderr, "[frag mini multi dp] node_i < 0, BUG.\n"); exit(0); }
                    line[0][node_i] = _left;
                    --node_i;
                    _right = _left;
                    _left = f_node[_right.x][_right.y].from;
                }
                if (node_i >= 0) { fprintf(stderr, "[frag mini multi dp] node_i >= 0, BUG.\n"); exit(0); }
                if (_right.x != left.x || _right.y != left.y)
                    line_end[0] = 0;
                else break;
            }
        }
    }
    //backtrack
    //Uni-Best  or  All-Best
    l_i = 1;
    for (i = 0; i < candi_n; ++i)
    {
        if (i > 0)	//check for the backtrack_flag
        {
            _right = candi[i];
            while (_right.x != head.x)
            {
                if (f_node[_right.x][0].backtrack_flag == 2)
                    goto NextLoop;
                _left = f_node[_right.x][_right.y].from;
                _right = _left;
            }
        }
        _right = candi[i];
        node_i = f_node[_right.x][_right.y].node_n - 1;
        line_end[l_i] = node_i+1;
        while (_right.x != head.x)
        {
            f_node[_right.x][0].backtrack_flag = 2;		//XXX
            if (node_i < 0) { fprintf(stderr, "[frag mini multi dp] node_i < 0, BUG.\n"); exit(0); }
            line[l_i][node_i] = _right;
            --node_i;
            _left = f_node[_right.x][_right.y].from;
            _right = _left;
        }
        if (node_i >= 0) { fprintf(stderr, "[frag mini multi dp] node_i >= 0, BUG.\n"); exit(0); }
        l_i++;
NextLoop:;
    }
    /*for (i = 1; i < l_i; ++i)
      {
      fprintf(stdout, "mini line# %d:\n", i);
      for (j = 0; j < line_end[i]; ++j)
      fprintf(stdout, "(%d,%d)\t", line[i][j].x, line[i][j].y);
      printf("\n");
      }*/
    int con_flag1, con_flag2, max_l;
    if (_head && _tail)	// choose the line that matches with anchor(s)
    {
        max_score = 0;
        max_l = 0;
        //there might be some overlaps????, might not be
        for (i = 1; i < l_i; ++i)
        {
            if (max_score > f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score)
                break;
            get_abs_dis(a_msg, left.x, left.y, line[i][0].x, line[i][0].y, &con_flag1, seed_len);	
            get_abs_dis(a_msg, right.x, right.y, line[i][line_end[i]-1].x, line[i][line_end[i]-1].y, &con_flag2, seed_len);
            if (con_flag1 == F_UNCONNECT || con_flag1 == F_CHR_DIF || con_flag2 == F_UNCONNECT || con_flag2 == F_CHR_DIF)
                continue;
            else if (f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score + F_SV_PEN + ((con_flag1==F_MATCH||con_flag1==F_MISMATCH)?F_MATCH_SCORE:-F_SV_PEN)
                    + 1 + ((con_flag2==F_MATCH||con_flag2==F_MISMATCH)?F_MATCH_SCORE:-F_SV_PEN) > max_score)
            {
                max_score = f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score + F_SV_PEN + ((con_flag1==F_MATCH||con_flag1==F_MISMATCH)?F_MATCH_SCORE:-F_SV_PEN)
                    + 1 + ((con_flag2==F_MATCH||con_flag2==F_MISMATCH)?F_MATCH_SCORE:F_SV_PEN);
                max_l = i;
            }
        }
        //move the match line to 'line[0]'
        if (max_l != 0)
        {
            for (i = 0; i < line_end[max_l]; ++i)
            {
                line[0][i] = line[max_l][i];
            }
            line_end[0] = line_end[max_l];
            line_end[max_l] = 0;
        }
        else line_end[0] = 0;
    }
    else if (_head)
    {
        max_score = 0;
        max_l = 0;
        for (i = 1; i < l_i; ++i)
        {
            if (max_score > f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score)
                break;
            get_abs_dis(a_msg, left.x, left.y, line[i][0].x, line[i][0].y, &con_flag1, seed_len);	
            if (con_flag1 == F_UNCONNECT || con_flag1 == F_CHR_DIF)
                continue;
            else if (f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score + F_SV_PEN +  ((con_flag1==F_MATCH||con_flag1==F_MISMATCH)?F_MATCH_SCORE:-F_SV_PEN) > max_score)
            {
                max_score = f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score + F_SV_PEN + ((con_flag1==F_MATCH||con_flag1==F_MISMATCH)?F_MATCH_SCORE:-F_SV_PEN);
                max_l = i;
            }
        }
        if (max_l != 0)
        {
            for (i = 0; i < line_end[max_l]; ++i)
            {
                line[0][i] = line[max_l][i];
            }
            line_end[0] = line_end[max_l];
            line_end[max_l] = 0;
        }
        else line_end[0] = 0;
    }
    else if (_tail)
    {
        max_score = 0;
        max_l = 0;
        for (i = 1; i < l_i; ++i)
        {
            if (max_score > f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score)
                break;
            get_abs_dis(a_msg, right.x, right.y, line[i][line_end[i]-1].x, line[i][line_end[i]-1].y, &con_flag2, seed_len);
            if (con_flag2 == F_UNCONNECT || con_flag2 == F_CHR_DIF)
                continue;
            else if (f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score + 1 + ((con_flag2==F_MATCH||con_flag2==F_MISMATCH)?F_MATCH_SCORE:-F_SV_PEN) > max_score)
            {
                max_score = f_node[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].score + 1 + ((con_flag2==F_MATCH||con_flag2==F_MISMATCH)?F_MATCH_SCORE:-F_SV_PEN);
                max_l = i;
            }
        }
        if (max_l != 0)
        {
            for (i = 0; i < line_end[max_l]; ++i)
            {
                line[0][i] = line[max_l][i];
            }
            line_end[0] = line_end[max_l];
            line_end[max_l] = 0;
        }
        else line_end[0] = 0;
    }
    return l_i;
}

/*int frag_mini_dp_multi_line(frag_dp_node **f_node, aln_msg *a_msg, int seed_len, line_node left, line_node right, line_node **line, int *line_end, int _head, int _tail, int n_seed)
{
    line_node head = ((_head)?left:(line_node){-1,0});
    int i, j, mini_dp_flag = MULTI_FLAG;
    //dp init
    for (i = left.x + 1; i < right.x; ++i)
    {
        //mini_dp XXX?
        //frag_dp_init(f_node, a_msg, i, left, seed_len, mini_dp_flag);
        frag_dp_multi_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);
    }
    //dp update
    for (i = left.x + 2; i < right.x; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            if (f_node[i][j].dp_flag == mini_dp_flag)
                frag_dp_multi_update(f_node, a_msg, i, j, left.x+1, seed_len, mini_dp_flag);
        }
    }
    //find backtrack start node
    int max_score = 0, max_n = 0;
    line_node max_node = head;//left;
    //UNLIMITED 
    if (_tail == 0)
    {
        for (i = right.x - 1; i > left.x; --i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                if (f_node[i][j].dp_flag == mini_dp_flag && f_node[i][j].score > max_score)
                {
                    max_score = f_node[i][j].score;
                    max_node = (line_node){i, j};
                    max_n = f_node[i][j].node_n;
                }
            }
        }
    }
    else
    {
        f_node[right.x][right.y].from = head;//left;
        f_node[right.x][right.y].score = 0;
        f_node[right.x][right.y].dp_flag = mini_dp_flag;
        f_node[right.x][right.y].node_n = 0;
        frag_dp_multi_update(f_node, a_msg, right.x, right.y, left.x+1, seed_len, mini_dp_flag);
        max_node = f_node[right.x][right.y].from;
        max_n = f_node[right.x][right.y].node_n;
    }
    //backtrack
    int l_i = 0;
    line_node _right = max_node, _left;
    int node_i = max_n-1;
    line_end[0] = node_i+1;
    //while (_right.x != left.x)
    while (_right.x != head.x)
    {
        if (node_i < 0) { fprintf(stderr, "[frag mini dp] node_i BUG 1.\n"); exit(0); }
        line[0][node_i--] = _right;
        _left = f_node[_right.x][_right.y].from;
        _right = _left;
    }
    if (node_i >= 0) { fprintf(stderr, "[frag mini dp] node_i BUG 2.\n"); exit(0); }
    ++l_i;

    //new lines
    if (right.x - left.x - 1 - max_n >= 2)
    {
        //dp init
        if (_head)
        {
            head.x = -1; head.y = 0;
            for (i = left.x + 1; i < right.x; ++i)
                frag_dp_multi_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);
        }
        for (i = 0; i < line_end[0]; ++i)
        {
            f_node[line[0][i].x][line[0][i].y].dp_flag = 0-mini_dp_flag;
            f_node[line[0][i].x][0].backtrack_flag = 2;	//indicate that this seed has been located already 
        }
        //dp update
        for (i = left.x + 2; i < right.x; ++i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                if (f_node[i][j].dp_flag == mini_dp_flag)
                    frag_dp_multi_update(f_node, a_msg, i, j, left.x+1, seed_len, mini_dp_flag);
            }
        }
        int k, l, candi_n=0, flag;
        line_node candi[n_seed];
        for (i = right.x-1; i > left.x; --i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                //candi-nodes
                if (f_node[i][j].backtrack_flag == 1 && f_node[i][j].dp_flag == mini_dp_flag && f_node[i][j].score > 0)
                {
                    //sort candi-nodes in 'candi', big->small
                    for (k = 0; k < candi_n; ++k)
                    {
                        if (f_node[i][j].score > f_node[line[0][k].x][line[0][k].y].score)
                            break;
                    }
                    for (l = candi_n; l > k; --l)
                    {
                        candi[l] = candi[l-1];
                    }
                    candi[k].x = i; candi[k].y = j;
                    candi_n++;
                }
            }
        }
        //XXX terminate early
        for (i = 0; i < candi_n; ++i)
        {
            //check for the backtrack_flag
            _right = candi[i];
            flag = 0;
            while (_right.x != head.x)
            {
                if (f_node[_right.x][0].backtrack_flag != 2)
                {
                    flag = 1;
                    break;
                }
                _left = f_node[_right.x][_right.y].from;
                _right = _left;
            }
            if (flag)
            {
                _right = candi[i];
                node_i = f_node[_right.x][_right.y].node_n - 1;
                line_end[l_i] = node_i+1;
                while (_right.x != head.x)
                {
                    f_node[_right.x][0].backtrack_flag = 2;		//XXX
                    if (node_i < 0) { fprintf(stderr, "[frag mini multi dp] node_i < 0, BUG.\n"); exit(0); }
                    line[l_i][node_i] = _right;
                    --node_i;
                    _left = f_node[_right.x][_right.y].from;
                    _right = _left;
                }
                if (node_i >= 0) { fprintf(stderr, "[frag mini multi dp] node_i >= 0, BUG.\n"); exit(0); }
                l_i++;
            }
        }
    }

    for (i = left.x+1; i < right.x; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            fprintf(stdout, "node:(%d %d)\t%d %d %d\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, a_msg[i].at[j].offset, f_node[i][j].from.x, f_node[i][j].from.y, f_node[i][j].score, "MXIDCRUSE"[f_node[i][j].match_flag], f_node[i][j].dp_flag, f_node[i][j].node_n);
        }
    }
    return l_i;
}*/

//for the uni-best line
/*int frag_mini_dp_line(frag_dp_node **f_node, 
        aln_msg *a_msg, 
        int seed_len, 
        line_node left, line_node right, 
        line_node *line, 
        int _head, int _tail, 
        int n_seed)
{
    line_node head = ((_head)?left:(line_node){-1,0});
    int i, j, mini_dp_flag = MULTI_FLAG;
    //dp init
    for (i = left.x + 1; i < right.x; ++i)
    {
        //mini_dp XXX?
        //frag_dp_init(f_node, a_msg, i, left, seed_len, mini_dp_flag);
        frag_dp_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);
    }
    //dp update
    for (i = left.x + 2; i < right.x; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            if (f_node[i][j].dp_flag == mini_dp_flag)
                frag_dp_update(f_node, a_msg, i, j, left.x+1, seed_len, mini_dp_flag);
        }
    }
    //find backtrack start node
    int max_score = 0, max_dis = 0, max_n = 0;
    line_node max_node = head;//left;
    //UNLIMITED 
    if (_tail == 0)
    {
        for (i = right.x - 1; i > left.x; --i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                if (f_node[i][j].dp_flag == mini_dp_flag)
                {
                    if (f_node[i][j].score > max_score)
                    {
                        max_score = f_node[i][j].score;
                        max_dis = f_node[i][j].dis_pen;
                        max_node = (line_node){i, j};
                        max_n = f_node[i][j].node_n;
                    }
                    else if (f_node[i][j].score == max_score)
                    {
                        if (f_node[i][j].dis_pen < max_dis)
                        {
                            max_score = f_node[i][j].score;
                            max_dis = f_node[i][j].dis_pen;
                            max_node = (line_node){i, j};
                            max_n = f_node[i][j].node_n;
                        }
                    }
                }
            }
        }
    }
    else
    {
        f_node[right.x][right.y].from = head;//left;
        f_node[right.x][right.y].score = 0;
        f_node[right.x][right.y].dp_flag = mini_dp_flag;
        f_node[right.x][right.y].node_n = 0;
        frag_dp_update(f_node, a_msg, right.x, right.y, left.x+1, seed_len, mini_dp_flag);
        max_node = f_node[right.x][right.y].from;
        max_n = f_node[right.x][right.y].node_n;
    }
    //backtrack
    line_node _right = max_node, _left;
    int node_i = max_n-1;
    //while (_right.x != left.x)
    while (_right.x != head.x)
    {
        if (node_i < 0) { fprintf(stderr, "[frag mini dp] node_i BUG 1.\n"); exit(0); }
        line[node_i--] = _right;
        _left = f_node[_right.x][_right.y].from;
        _right = _left;
    }
    if (node_i >= 0) { fprintf(stderr, "[frag mini dp] node_i BUG 2.\n"); exit(0); }

    for (i = left.x+1; i < right.x; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            fprintf(stdout, "node:(%d %d)\t%d %d %d\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, a_msg[i].at[j].offset, f_node[i][j].from.x, f_node[i][j].from.y, f_node[i][j].score, "MXIDCRUSE"[f_node[i][j].match_flag], f_node[i][j].dp_flag, f_node[i][j].node_n);
        }
    }
    return max_n;
}*/

//for best coverage and connect
int frag_mini_dp_line(frag_dp_node **f_node, 
                      aln_msg *a_msg, 
                      int seed_len, 
                      line_node left, line_node right, 
                      line_node *line, 
                      int _head, int _tail, 
                      int n_seed)
{
    line_node head = ((_head)?left:(line_node){-1,0});
    int i, j, mini_dp_flag = MULTI_FLAG;
    //dp init
    for (i = left.x + 1; i < right.x; ++i)
    {
        //mini_dp XXX?
        //frag_dp_init(f_node, a_msg, i, left, seed_len, mini_dp_flag);
        frag_dp_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);
    }
    //dp update
    for (i = left.x + 2; i < right.x; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            if (f_node[i][j].dp_flag == mini_dp_flag)
                frag_dp_update(f_node, a_msg, i, j, left.x+1, seed_len, mini_dp_flag);
        }
    }
    //find backtrack start node
    int max_score = 0, max_dis = 0, max_n = 0;
    line_node max_node = head;//left;
    //UNLIMITED 
    if (_tail == 0)
    {
        for (i = right.x - 1; i > left.x; --i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                if (f_node[i][j].dp_flag == mini_dp_flag)
                {
                    if (f_node[i][j].score > max_score)
                    {
                        max_score = f_node[i][j].score;
                        max_dis = f_node[i][j].dis_pen;
                        max_node = (line_node){i, j};
                        max_n = f_node[i][j].node_n;
                    }
                    else if (f_node[i][j].score == max_score)
                    {
                        if (f_node[i][j].dis_pen < max_dis)
                        {
                            max_score = f_node[i][j].score;
                            max_dis = f_node[i][j].dis_pen;
                            max_node = (line_node){i, j};
                            max_n = f_node[i][j].node_n;
                        }
                    }
                }
            }
        }
    }
    else
    {
        f_node[right.x][right.y].from = head;//left;
        f_node[right.x][right.y].score = 0;
        f_node[right.x][right.y].dp_flag = mini_dp_flag;
        f_node[right.x][right.y].node_n = 0;
        //frag_dp_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);
        f_node[right.x][right.y].node_n = 0;
        frag_dp_update(f_node, a_msg, right.x, right.y, left.x+1, seed_len, mini_dp_flag);
        max_node = f_node[right.x][right.y].from;
        max_n = f_node[right.x][right.y].node_n;
    }
    //backtrack
    line_node _right = max_node, _left;
    int node_i = max_n-1;
    //while (_right.x != left.x)
    while (_right.x != head.x)
    {
        if (node_i < 0) { fprintf(stderr, "[frag mini dp] node_i BUG 1.\n"); exit(0); }
        line[node_i--] = _right;
        _left = f_node[_right.x][_right.y].from;
        _right = _left;
    }
    if (node_i >= 0) { fprintf(stderr, "[frag mini dp] node_i BUG 2.\n"); exit(0); }

    /*for (i = left.x+1; i < right.x; ++i)
    {
        for (j = 0; j < a_msg[i].n_aln; ++j)
        {
            fprintf(stdout, "node:(%d %d)\t%d %d %d\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, a_msg[i].at[j].offset, f_node[i][j].from.x, f_node[i][j].from.y, f_node[i][j].score, "MXIDCRUSE"[f_node[i][j].match_flag], f_node[i][j].dp_flag, f_node[i][j].node_n);
        }
    }*/
    return max_n;
}

int frag_dp_line(aln_msg *a_msg, 
                 int n_seed, int seed_len, 
                 line_node **line, int *line_end, 
                 frag_dp_node ***f_node, 
                 line_node **_line, int *_line_end) //_headi=0, _tail=0
{
    int i, j, k;
    int min_len = 1, min_exist=0;

    //dp init
    {
        for (i = 0; i < n_seed; ++i)
        {
            if (a_msg[i].n_aln <= min_len)
            {
                frag_dp_init(*f_node, a_msg, i, (line_node){-1,0}, seed_len, MIN_FLAG);
                min_exist = 1;
            }
            else
                frag_dp_init(*f_node, a_msg, i, (line_node){-1,0}, seed_len, MULTI_FLAG);
        }
    }
    //dp update and backtrack
    {
        if (min_exist)
        {
            //min extend XXX
            /*for (i = 0; i < n_seed; ++i)
              {
              if (a_msg[i].n_aln <= min_len)
              {
              for (j = 0; j < a_msg[i].n_aln; ++j)
              frag_min_extend(f_node, i, j, n_seed, min_len, MIN_FLAG, seed_len);
              }
              }*/
            //min update
            for (i = 1; i < n_seed; ++i)
            {
                for (j = 0; j < a_msg[i].n_aln; ++j)
                {
                    if ((*f_node)[i][j].dp_flag == MIN_FLAG)
                        frag_dp_update(*f_node, a_msg, i, j, 0/*update start pos*/, seed_len, MIN_FLAG);
                }
            }
            //find start node of backtrack
            int max_score=0, max_dis = 0, l_i = 0;
            line_node max_node = (line_node){-1,0};
            int node_i=0, mini_len, new_l=1;
            line_node last_n; int multi_l, max_multi = n_seed > 10 ? 10 : n_seed;
            line_node right, left;
            for (i = n_seed-1; i >= 0; --i)
            {
                for (j = 0; j < a_msg[i].n_aln; ++j)
                {
                    if ((*f_node)[i][j].backtrack_flag == 1 && (*f_node)[i][j].dp_flag == MIN_FLAG && (*f_node)[i][j].score > 0)
                    {
                        if ((*f_node)[i][j].score > max_score)
                        {
                            max_score = (*f_node)[i][j].score;
                            max_node = (line_node){i,j};
                            max_dis = (*f_node)[i][j].dis_pen;
                        }
                        else if ((*f_node)[i][j].score == max_score)
                        {
                            if ((*f_node)[i][j].dis_pen < max_dis)
                            {
                                max_score = (*f_node)[i][j].score;
                                max_node = (line_node){i,j};
                                max_dis = (*f_node)[i][j].dis_pen;
                            }
                        }
                    }
                }
            }
            //backtrack
            //ONLY one backtrack node.
            if (max_node.x < n_seed - 1)
            {
                //mini-dp with anchors
                mini_len = frag_mini_dp_line(*f_node, a_msg, seed_len, max_node, (line_node){n_seed, 0}, _line[0], 1, 0, n_seed);
                for (k = mini_len - 1; k >= 0; --k)
                    line[l_i][node_i++] = _line[0][k];	//_line[0][k] is the last node
                line[l_i][node_i] = max_node;	//twice write, for the last multi-dp-line	//XXX
                //add multi-line
                last_n = (line_node){n_seed, 0};
                for (k = mini_len; k >= 0; --k)
                {
                    if (last_n.x - line[l_i][node_i-k].x > 2)
                    {
                        //add multi-line
                        multi_l = frag_mini_dp_multi_line(*f_node, a_msg, seed_len, line[l_i][node_i-k], last_n, _line, _line_end, max_multi, 0, 0, n_seed);
                        for (i = 1; i < multi_l; ++i)
                        {
                            if (_line_end[i] > 0)
                            {
                                for (j = 0; j < _line_end[i]; ++j)
                                    line[l_i+new_l][j] = _line[i][j];
                                line_end[l_i+new_l] = _line_end[i];
                                line[l_i+new_l][line_end[l_i+new_l]] = line[l_i][node_i-k];
                                line[l_i+new_l][line_end[l_i+new_l]+1] = last_n;
                                line[l_i+new_l][line_end[l_i+new_l]+2] = (line_node){0, 0};	//need trigger
                                (*f_node)[last_n.x][last_n.y].trigger[(*f_node)[last_n.x][last_n.y].next_trigger_n++] = l_i+new_l;
                                ++new_l;
                            }
                        }
                    }
                    last_n = line[l_i][node_i-k];
                }
                /*mini_l = frag_mini_dp_multi_line(f_node, a_msg, seed_len, max_node, (line_node){n_seed, 0}, _line, _line_end, 1, 0, n_seed);
                //might be cut off or not
                for (k = _line_end[0]-1; k >= 0; --k)
                line[l_i][node_++] = _line[0][k];
                for (i = 1; i < mini_l; ++i)
                {
                if (_line_end[i] > 0)
                {
                for (j = 0; j < _line_end[i]; ++j)
                {
                line[l_i+new_l][j] = _line[i][j];
                }
                line_end[l_i+new_l] = _line_end[i];
                line[l_i+new_l][line_end[l_i+new_l]] = max_node;
                line[l_i+new_l][line_end[l_i+new_l]+1] = (line_node){n_seed, 0};
                ++new_l;
                }
                }*/
            }
            right = max_node;
            while (right.x != -1)
            {
                line[l_i][node_i++] = right;
                left = (*f_node)[right.x][right.y].from;

                if (left.x < right.x - 1 )//f_node[right.x][right.y].match_flag != F_MATCH)	//XXX if left.x < right.x-1, match_flag of right node couldn't be match?
                {
                    //mini-dp with anchors
                    mini_len = frag_mini_dp_line(*f_node, a_msg, seed_len, left, right, _line[0], 1, 1, n_seed);
                    for (k = mini_len-1; k >= 0; --k)
                        line[l_i][node_i++] = _line[0][k];
                    line[l_i][node_i] = left;	//twice write, for the last multi-dp-line	//XXX
                    //add multi-line
                    last_n = right;
                    for (k = mini_len; k >= 0; --k)
                    {
                        if (last_n.x - line[l_i][node_i-k].x > 2)
                        {
                            //add multi-line
                            multi_l = frag_mini_dp_multi_line(*f_node, a_msg, seed_len, line[l_i][node_i-k], last_n, _line, _line_end, max_multi, 0, 0, n_seed);
                            for (i = 1; i < multi_l; ++i)
                            {
                                if (_line_end[i] > 0)
                                {
                                    for (j = 0; j < _line_end[i]; ++j)
                                        line[l_i+new_l][j] = _line[i][j];
                                    line_end[l_i+new_l] = _line_end[i];
                                    line[l_i+new_l][line_end[l_i+new_l]] = line[l_i][node_i-k];
                                    line[l_i+new_l][line_end[l_i+new_l]+1] = last_n;
                                    line[l_i+new_l][line_end[l_i+new_l]+2] = (line_node){0, 0};	//need trigger
                                    (*f_node)[last_n.x][last_n.y].trigger[(*f_node)[last_n.x][last_n.y].next_trigger_n + (*f_node)[last_n.x][last_n.y].pre_trigger_n] = l_i+new_l;
                                    ++(*f_node)[last_n.x][last_n.y].pre_trigger_n;
                                    ++new_l;
                                }
                            }
                        }
                        last_n = line[l_i][node_i-k];
                    }
                    /*if (left.x == -1)
                      mini_l = frag_mini_dp_multi_line(f_node, a_msg, seed_len, left, right, _line, _line_end, 0, 1, n_seed);
                      else
                      mini_l = frag_mini_dp_multi_line(f_node, a_msg, seed_len, left, right, _line, _line_end, 1, 1, n_seed);
                    //fprintf(stdout, "left: %d, right: %d\n", left.x, right.x);
                    //might be cut off or not
                    for (k = _line_end[0]-1; k >= 0; --k)
                    line[l_i][node_i++] = _line[0][k];
                    for (i = 1; i < mini_l; ++i)
                    {
                    if (_line_end[i] > 0)
                    {
                    for (j = 0; j < _line_end[i]; ++j)
                    {
                    line[l_i+new_l][j] = _line[i][j];
                    }
                    line_end[l_i+new_l] = _line_end[i];
                    line[l_i+new_l][line_end[l_i+new_l]] = left;
                    line[l_i+new_l][line_end[l_i+new_l]+1] = right;
                    ++new_l;
                    }
                    }*/
                }
                right = left;
            }
            //invert line
            line_node tmp;
            for (k = 0; k < node_i/2; ++k)
            {
                tmp =line[l_i][k]; line[l_i][k] = line[l_i][node_i-k-1]; line[l_i][node_i-k-1] = tmp;
            }

            //print
            /*for (i = 0; i < n_seed; ++i)
              {
              for (j = 0; j < a_msg[i].n_aln; ++j)
              {
              fprintf(stdout, "node:(%d %d)\t%d %d %d\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, a_msg[i].at[j].offset, f_node[i][j].from.x, f_node[i][j].from.y, f_node[i][j].score, "MXIDCRUSE"[f_node[i][j].match_flag], f_node[i][j].dp_flag, f_node[i][j].node_n);
              }
              }*/
            line_end[l_i] = node_i;
            line[l_i][line_end[l_i]] = (line_node){-1,0};
            line[l_i][line_end[l_i]+1] = (line_node){n_seed, 0};
            line[l_i][line_end[l_i]+2] = (line_node){1,1};	//trigger flag: (1,1) -> no needs for trigger
            //				(0,0) -> need trigger
            /*for (i = 0; i < new_l; ++i)
              {
              fprintf(stdout, "line #%d:\n", i+1);
              for (j = 0; j < line_end[i]+2; ++j)
              {
              fprintf(stdout, "(%d,%d)\t", line[i][j].x, line[i][j].y);
              }
              fprintf(stdout, "\n");
              }*/
            return l_i+new_l;
        }
        else	//whole-multi dp, should not go to here for current data.
        {
            fprintf(stderr, "[frag dp line] No min-seed.\n");
            return 0;
        }
    }
}

int frag_min_extend(frag_dp_node **f_node, aln_msg *a_msg,
                    int node_i, int aln_i,
                    int node_n, int aln_min,
                    int dp_flag,
                    int seed_len)
{
    int i, j, con_flag;
    int last_x = node_i, last_y = aln_i;

    // 0 -> node_i-1, node_i+1 -> node_n-1
    //from right to left 
    i = node_i - 1;
    while (i >= 0)
    {
        if (a_msg[i].n_aln > aln_min)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
            //dis = get_abs_dis(a_msg, from.x, from.y, seed_i, i, &con_flag, seed_len);
                get_abs_dis(a_msg, i, j, node_i, aln_i, &con_flag, seed_len);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH)
                {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        --i;
    }

    i = node_i + 1;
    while (i < node_n)
    {
        if (a_msg[i].n_aln > aln_min)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                get_abs_dis(a_msg, node_i, aln_i, i, j, &con_flag, seed_len);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH)
                {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        ++i;
    }
    return 0;
}

//Best Coverage and Connect
int frag_line_BCC(aln_msg *a_msg,
                  int n_seed, int seed_len,
                  line_node **line, int *line_end,
                  frag_dp_node ***f_node,
                  line_node **_line, int *_line_end)
{
    int i, j, k;
    // min_len XXX 
    int min_len = 1, min_exist=0, min_num = 0;

    //dp init
    {
        for (i = 0; i < n_seed; ++i)
        {
            if (a_msg[i].n_aln <= min_len)
            {
                frag_dp_init(*f_node, a_msg, i, (line_node){-1,0}, seed_len, MIN_FLAG);
                min_exist = 1;
                ++min_num;
            }
            else
                frag_dp_init(*f_node, a_msg, i, (line_node){-1,0}, seed_len, MULTI_FLAG);
        }
        //fraction XXX
        if (!min_exist || min_num * 10 < n_seed)
        {
            for (i = 0; i < n_seed; ++i)
            {
                //frag_dp_init(*f_node, a_msg, i, (line_node){-1,0}, seed_len, MIN_FLAG);
                for (j = 0; j < a_msg[i].n_aln; ++j) (*f_node)[i][j].dp_flag = MIN_FLAG;
            }
            min_len = PER_ALN_N; // XXX
            min_exist = 1;
        }
    }
    //dp update and backtrack
    {
        //min extend XXX, when min_len == PER_ALN_N: no need to extend
        for (i = 0; i < n_seed; ++i)
        {
            if (a_msg[i].n_aln <= min_len)
            {
                for (j = 0; j < a_msg[i].n_aln; ++j)
                    frag_min_extend(*f_node, a_msg, i, j, n_seed, min_len, MIN_FLAG, seed_len);
            }
        }
        //min update
        for (i = 1; i < n_seed; ++i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                if ((*f_node)[i][j].dp_flag == MIN_FLAG)
                    frag_dp_update(*f_node, a_msg, i, j, 0/*update start pos*/, seed_len, MIN_FLAG);
            }
        }
        //print
        /*printf("min-update:\n");
        for (i = 0; i < n_seed; ++i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
                //if ((*f_node)[i][j].dp_flag == MIN_FLAG)
                    fprintf(stdout, "node:(%d %d)\t%d\t%d %d %lld\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].read_id, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, (long long)a_msg[i].at[j].offset, (*f_node)[i][j].from.x, (*f_node)[i][j].from.y, (*f_node)[i][j].score, FRAG_CON_STR[(*f_node)[i][j].match_flag], (*f_node)[i][j].dp_flag, (*f_node)[i][j].node_n);
        }*/
        //find start node of backtrack
        int max_score=0, max_dis = 0, l_i = 0;
        line_node max_node = (line_node){-1,0};
        int node_i=0, mini_len, new_l=1;
        line_node last_n; int multi_l, max_multi = n_seed > 10 ? 10 : n_seed; //XXX
        line_node right, left;
        for (i = n_seed-1; i >= 0; --i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                if ((*f_node)[i][j].backtrack_flag == 1 && (*f_node)[i][j].dp_flag == MIN_FLAG && (*f_node)[i][j].score > 0)
                {
                    if ((*f_node)[i][j].score > max_score)
                    {
                        max_score = (*f_node)[i][j].score;
                        max_node = (line_node){i,j};
                        max_dis = (*f_node)[i][j].dis_pen;
                    }
                    else if ((*f_node)[i][j].score == max_score)
                    {
                        if ((*f_node)[i][j].dis_pen < max_dis)
                        {
                            max_score = (*f_node)[i][j].score;
                            max_node = (line_node){i,j};
                            max_dis = (*f_node)[i][j].dis_pen;
                        }
                    }
                }
            }
        }
        //backtrack
        if (max_node.x == -1) return 0;

        int tri_x, tri_y;
        //ONLY one backtrack node.
        if (max_node.x < n_seed - 1)
        {
            //mini-dp with anchors
            mini_len = frag_mini_dp_line(*f_node, a_msg, seed_len, max_node, (line_node){n_seed, 0}, _line[0], 1, 0, n_seed);

            if (mini_len > 0)
            {
                for (k = mini_len - 1; k >= 0; --k)
                    line[l_i][node_i++] = _line[0][k];	//_line[0][k] is the last node
            }


            line[l_i][node_i] = max_node;	//twice write, for the last multi-dp-line	//XXX
            //add trigger for multi-line
            last_n = (line_node){n_seed, 0};
            for (k = mini_len; k >= 0; --k)
            {
                if (last_n.x - line[l_i][node_i-k].x > 2)
                {
                    tri_x = line[l_i][node_i-k].x;
                    tri_y = line[l_i][node_i-k].y;
                    //add multi-line
                    multi_l = frag_mini_dp_multi_line(*f_node, a_msg, seed_len, line[l_i][node_i-k], last_n, _line, _line_end, max_multi, 0, 0, n_seed);
                    for (i = 1; i < multi_l; ++i)
                    {
                        if (_line_end[i] > 0)
                        {
                            for (j = 0; j < _line_end[i]; ++j)
                                line[l_i+new_l][j] = _line[i][j];
                            line_end[l_i+new_l] = _line_end[i];
                            line[l_i+new_l][line_end[l_i+new_l]] = line[l_i][node_i-k];
                            line[l_i+new_l][line_end[l_i+new_l]+1] = last_n;
                            line[l_i+new_l][line_end[l_i+new_l]+2] = (line_node){0, 0};	//need trigger
                            if ((*f_node)[tri_x][tri_y].next_trigger_n == (*f_node)[tri_x][tri_y].trigger_m)
                            {
                                (*f_node)[tri_x][tri_y].trigger_m <<= 1;
                                (*f_node)[tri_x][tri_y].trigger = (int*)realloc((*f_node)[tri_x][tri_y].trigger,  (*f_node)[tri_x][tri_y].trigger_m * sizeof(int));
                            }
                            (*f_node)[tri_x][tri_y].trigger[(*f_node)[tri_x][tri_y].next_trigger_n++] = l_i+new_l;
                            ++new_l;
                        }
                    }
                }
                last_n = line[l_i][node_i-k];
            }
        }
        right = max_node;
        while (right.x != -1)
        {
            line[l_i][node_i++] = right;
            left = (*f_node)[right.x][right.y].from;

            if (left.x < right.x - 1 )//f_node[right.x][right.y].match_flag != F_MATCH)	//XXX if left.x < right.x-1, match_flag of right node couldn't be match?
            {
                //mini-dp with anchors
                mini_len = frag_mini_dp_line(*f_node, a_msg, seed_len, left, right, _line[0], 1, 1, n_seed);
                for (k = mini_len-1; k >= 0; --k)
                    line[l_i][node_i++] = _line[0][k];
                line[l_i][node_i] = left;	//twice write, for the last multi-dp-line	//XXX
                // add trigger multi-len
                last_n = right;
                for (k = mini_len; k >= 0; --k)
                {
                    if (last_n.x - line[l_i][node_i-k].x > 2)
                    {
                        //add multi-line
                        multi_l = frag_mini_dp_multi_line(*f_node, a_msg, seed_len, line[l_i][node_i-k], last_n, _line, _line_end, max_multi, 0, 0, n_seed);
                        for (i = 1; i < multi_l; ++i)
                        {
                            if (_line_end[i] > 0)
                            {
                                for (j = 0; j < _line_end[i]; ++j)
                                    line[l_i+new_l][j] = _line[i][j];
                                line_end[l_i+new_l] = _line_end[i];
                                line[l_i+new_l][line_end[l_i+new_l]] = line[l_i][node_i-k];
                                line[l_i+new_l][line_end[l_i+new_l]+1] = last_n;
                                line[l_i+new_l][line_end[l_i+new_l]+2] = (line_node){0, 0};	//need trigger
                                if ((*f_node)[last_n.x][last_n.y].next_trigger_n + (*f_node)[last_n.x][last_n.y].pre_trigger_n == (*f_node)[last_n.x][last_n.y].trigger_m)
                                {
                                    (*f_node)[last_n.x][last_n.y].trigger_m <<= 1;
                                    (*f_node)[last_n.x][last_n.y].trigger = (int*)realloc((*f_node)[last_n.x][last_n.y].trigger,  (*f_node)[last_n.x][last_n.y].trigger_m * sizeof(int));
                                }
                                (*f_node)[last_n.x][last_n.y].trigger[(*f_node)[last_n.x][last_n.y].next_trigger_n + (*f_node)[last_n.x][last_n.y].pre_trigger_n] = l_i+new_l;
                                ++(*f_node)[last_n.x][last_n.y].pre_trigger_n;
                                ++new_l;
                            }
                        }
                    }
                    last_n = line[l_i][node_i-k];
                }
            }
            right = left;
        }
        //invert line
        line_node tmp;
        for (k = 0; k < node_i/2; ++k)
        {
            tmp =line[l_i][k]; line[l_i][k] = line[l_i][node_i-k-1]; line[l_i][node_i-k-1] = tmp;
        }

        //print
        /*printf("whole-update:\n");
        for (i = 0; i < n_seed; ++i)
        {
            for (j = 0; j < a_msg[i].n_aln; ++j)
            {
                fprintf(stdout, "node:(%d %d)\t%d\t%d %d %lld\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tBK-flag:%d\tnode_n:%d\n", i, j, a_msg[i].read_id, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, (long long)a_msg[i].at[j].offset, (*f_node)[i][j].from.x, (*f_node)[i][j].from.y, (*f_node)[i][j].score, FRAG_CON_STR[(*f_node)[i][j].match_flag], (*f_node)[i][j].dp_flag, (*f_node)[i][j].backtrack_flag, (*f_node)[i][j].node_n);
            }
        }*/
        line_end[l_i] = node_i;
        line[l_i][line_end[l_i]] = (line_node){-1,0};
        line[l_i][line_end[l_i]+1] = (line_node){n_seed, 0};
        line[l_i][line_end[l_i]+2] = (line_node){1,1};	//trigger flag: (1,1) -> no needs for trigger
                                                        //				(0,0) -> need trigger
        
        /*for (i = 0; i < new_l; ++i)
        {
            fprintf(stdout, "line #%d:\n", i+1);
            for (j = 0; j < line_end[i]+3; ++j)
            {
                fprintf(stdout, "(%d,%d)\t", line[i][j].x, line[i][j].y);
            }
            fprintf(stdout, "\n");
        }*/
        return l_i+new_l;
    }
}
//all valid aln 

int frag_mini_dp_path(aln_msg *a_msg, 
        int seed_len, 
        frag_msg *f_msg, 
        frag_dp_node **f_node, 
        line_node *line, 
        int m_len)
{
    int i;
    int frag_num = 0;
    int cur_x, cur_y, pre_x = line[m_len-1].x, pre_y = line[m_len-1].y;

    frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, f_msg, frag_num, seed_len);
    for (i = m_len-1; i > 0; --i)
    {
        cur_x = pre_x; cur_y = pre_y;
        pre_x = line[i-1].x; pre_y = line[i-1].y;

        if (f_node[cur_x][cur_y].match_flag != F_MATCH && f_node[cur_x][cur_y].match_flag != F_MISMATCH)
        {
            frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, f_msg, frag_num, seed_len);
            ++frag_num;
            frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, f_msg, frag_num, seed_len);
        }
        else
            frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, f_msg, frag_num, seed_len);
    }
    cur_x = line[0].x; cur_y = line[0].y;
    frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, f_msg, frag_num, seed_len);
    return 1;
}

//@para:
//	a_msg:	struct of seeds' aln result
//	n_seed:	whole number of seeds that have at least 1 aln-res, and less than 100
//	f_msg:	be used to store frag-msg
//	line_n:	number of lines
//	line_m:	max number of lines allowed in mem
int _frag_dp_path(aln_msg *a_msg, 
        int n_seed, int seed_len, 
        frag_msg **f_msg, 
        int *line_n/*return, line_num*/, int *line_tri, int *line_m/*return, line mem size*/, 
        line_node **line, int *line_end, 
        frag_dp_node ***f_node, 
        line_node **_line, int *_line_end,
        int aln_type)
{
    int i, l;
    //DP
    int l_n;
    if (aln_type == 1) // best coverage and connect
        l_n = frag_line_BCC(a_msg, n_seed, seed_len, line, line_end, f_node, _line, _line_end);// for best coverage and connect case: l_n equals 1.
    else l_n = 0; // all valid
    if (l_n == 0) return 0;

    /*for (i = 0; i < l_n; ++i)
      {
      printf("%d:\t", i+1);
      for (l = 0; l < line_end[i]; ++l)
      printf("(%d, %d)\t", a_msg[line[i][l].x].read_id, line[i][l].y);
      printf("\n");
      }*/

    int frag_num;
    int cur_x , cur_y , pre_x, pre_y;

    //for new frags
    //  cur_line: cur line index, cur_num:  whole lines, *line_num: mem of *f_msg
    int cur_line=0;//, cur_num=1, _m_len;
    int left_bound, right_bound;

    // calcu whole number of line
    (*line_n) = 1;
    for (i = line_end[0]-1; i > 0; --i)
    {
        if ((*f_node)[line[0][i].x][line[0][i].y].match_flag == F_CHR_DIF || (*f_node)[line[0][i].x][line[0][i].y].match_flag == F_UNCONNECT)
            ++(*line_n);
    }
    if (*line_n > (*line_m))
    {
        (*f_msg) = (frag_msg*)realloc(*f_msg, *line_n * sizeof(frag_msg));
        for (i = (*line_m); i < *line_n; ++i)
        {
            frag_init_msg((*f_msg)+i, (*f_msg)->frag_max); 
            frag_copy_msg(*f_msg, (*f_msg)+i);
        }
        if ((*f_msg) == NULL) { fprintf(stderr, "[frag_dp_path] Not enough memory.(line_m: %d)\n", *line_n); exit(0); }
        (*line_m)= *line_n;
        fprintf(stderr, "line-num: %d\t", *line_n);
    }

    for (l = 0; l < l_n; ++l)
    {
        //line_tri[l] = line[l][line_end[l]+2].x;
        frag_num=0;
        pre_x = line[l][line_end[l]-1].x;
        pre_y = line[l][line_end[l]-1].y;
        //fprintf(stdout, "left: %d, right: %d\n", line[l][line_end[l]].x, line[l][line_end[l]+1].x);
        //right bound
        right_bound = ((line[l][line_end[l]+1].x == n_seed) ? ((*f_msg)[0].seed_all+1) : (a_msg[line[l][line_end[l]+1].x].read_id));
        //first end
        frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
        ((*f_msg)+cur_line)->frag_right_bound = right_bound;
        line_tri[cur_line] = 1;
        for (i = line_end[l]-1; i > 0; --i)
        {
            cur_x = pre_x; cur_y = pre_y;
            pre_x = line[l][i-1].x; pre_y = line[l][i-1].y;
            //INS
            if ((*f_node)[cur_x][cur_y].match_flag == F_INSERT)
            {
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
                ++frag_num;
                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
            }
            //DEL
            else if ((*f_node)[cur_x][cur_y].match_flag == F_DELETE)
            {
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
                ++frag_num;
                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
            }
            //MIS	XXX
            else if ((*f_node)[cur_x][cur_y].match_flag == F_MISMATCH || (*f_node)[cur_x][cur_y].match_flag == F_LONG_MISMATCH)
            {
                //XXX take mis-match as NEW-FLAG case
                //frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
                ++frag_num;
                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
            }
            else if ((*f_node)[cur_x][cur_y].match_flag == F_MATCH)	//: MATCH
            {
                frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
            }
            else if ((*f_node)[cur_x][cur_y].match_flag == F_CHR_DIF) // : INV
            {
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
                left_bound = a_msg[pre_x].read_id;
                ((*f_msg)+cur_line)->frag_left_bound = left_bound;

                right_bound = a_msg[cur_x].read_id;
                ++cur_line;
                if (cur_line >= (*line_m)) fprintf(stderr, "[frag dp path] line number error.\n");
                frag_num = 0;
                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
                ((*f_msg)+cur_line)->frag_right_bound = right_bound;
                line_tri[cur_line] = 1;
            }
            else if ((*f_node)[cur_x][cur_y].match_flag == F_UNCONNECT) // : Alu ...
            {
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
                left_bound = a_msg[pre_x].read_id;
                ((*f_msg)+cur_line)->frag_left_bound = left_bound;

                right_bound = a_msg[cur_x].read_id;
                ++cur_line;
                if (cur_line >= (*line_m)) fprintf(stderr, "[frag dp path] line number error.\n");
                frag_num = 0;
                frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
                ((*f_msg)+cur_line)->frag_right_bound = right_bound;
                line_tri[cur_line] = 1;
            }
            else {fprintf(stderr, "[frag dp path] Error: Unknown flag, \"%d\"", (*f_node)[cur_x][cur_y].match_flag); exit(0);}
        }
        //last start
        cur_x = line[l][0].x; cur_y = line[l][0].y;
        frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
        left_bound = ((line[l][line_end[l]].x == -1) ? 0 : (a_msg[line[l][line_end[l]].x].read_id));
        ((*f_msg)+cur_line)->frag_left_bound = left_bound;
    }

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

int frag_dp_path(aln_msg *a_msg, 
				 int n_seed, int seed_len, 
				 frag_msg **f_msg, 
				 int *line_n, int *line_tri, int *line_m, 
				 line_node **line, int *line_end, 
				 frag_dp_node ***f_node, 
				 line_node **_line, int *_line_end,
                 int aln_type)
{
	int i, l;
	//DP
    int l_n;
    if (aln_type == 1) // best coverage 
        l_n = frag_line_BCC(a_msg, n_seed, seed_len, line, line_end, f_node, _line, _line_end);
    //else
	if (l_n == 0) return 0;

    /*for (i = 0; i < l_n; ++i)
    {
        printf("%d:\t", i+1);
        for (l = 0; l < line_end[i]; ++l)
            printf("(%d, %d)\t", a_msg[line[i][l].x].read_id, line[i][l].y);
        printf("\n");
    }*/

	int frag_num;
	int cur_x , cur_y , pre_x, pre_y;

	//for new frags
	//  cur_line: cur line index, cur_num:  whole lines, *line_num: mem of *f_msg
	//int cur_line=0, cur_num=1, _m_len;
	int left_bound, right_bound;

    (*line_n) = l_n;
    if (l_n > (*line_m))
    {
        (*f_msg) = (frag_msg*)realloc(*f_msg, l_n * sizeof(frag_msg));
        for (i = (*line_m); i < l_n; ++i)
        {
            frag_init_msg((*f_msg)+i, (*f_msg)->frag_max); 
            frag_copy_msg(*f_msg, (*f_msg)+i);
        }
        if ((*f_msg) == NULL) { fprintf(stderr, "[frag_dp_path] Not enough memory.(line_m: %d)\n", l_n); exit(0); }
        (*line_m)= l_n;
        fprintf(stderr, "line-num: %d\t", l_n);
    }

	for (l = 0; l < l_n; ++l)
	{
        line_tri[l] = line[l][line_end[l]+2].x;
        frag_num=0;
        pre_x = line[l][line_end[l]-1].x;
		pre_y = line[l][line_end[l]-1].y;
        //fprintf(stdout, "left: %d, right: %d\n", line[l][line_end[l]].x, line[l][line_end[l]+1].x);
		//right bound
        right_bound = ((line[l][line_end[l]+1].x == n_seed) ? ((*f_msg)[0].seed_all+1) : (a_msg[line[l][line_end[l]+1].x].read_id));
		//right_bound = (*f_msg)[0].seed_all+1;
		//MIS-MATCH
		/*if (abs(line[l][line_end[l]+1].x-pre_x) > 1)
		{
			_m_len = 0;
			_m_len = frag_mini_dp_multi_line((*f_node), a_msg, seed_len, (line_node){pre_x, pre_y}, (line_node){n_seed, 0}, _line, _line_end, 0, 0, n_seed);
			if (_m_len != 0)
			for (i = 1; i < _m_len; ++i)
			{
				if (cur_num == (*line_m))
				{
					(*f_msg) = (frag_msg*)realloc(*f_msg, ((*line_m)+1) * sizeof(frag_msg));
					frag_init_msg((*f_msg)+(*line_m), (*f_msg)->frag_max); 
					if ((*f_msg) == NULL) { fprintf(stderr, "[frag_dp_path] Not enough memory.(line_m: %d)\n", (*line_m)+1); exit(0); }
					(*line_m)++;
				}
				frag_copy_msg(*f_msg, (*f_msg)+cur_num);
				frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, *f_node, _line[i], _line_end[i]);
				((*f_msg)+cur_num)->frag_left_bound = a_msg[pre_x].read_id;
				((*f_msg)+cur_num)->frag_right_bound = (*f_msg)[0].seed_all+1;
				cur_num++;
				//right_bound = a_msg[_line[i][0].x].read_id;	//XXX
			}
		}*/
		//first end
		frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+l, frag_num, seed_len);
		if ((*f_node)[pre_x][pre_y].next_trigger_n > 0 || (*f_node)[pre_x][pre_y].pre_trigger_n > 0)	//set trigger
            frag_trigger_set((*f_node)[pre_x][pre_y], (*f_msg)+l, frag_num);
        ((*f_msg)+l)->frag_right_bound = right_bound;
		for (i = line_end[l]-1; i > 0; --i)
		{
			cur_x = pre_x; cur_y = pre_y;
			pre_x = line[l][i-1].x; pre_y = line[l][i-1].y;
			//INS
			if ((*f_node)[cur_x][cur_y].match_flag == F_INSERT)
			{
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+l, frag_num, seed_len);
				if ((*f_node)[cur_x][cur_y].next_trigger_n > 0 || (*f_node)[cur_x][cur_y].pre_trigger_n > 0)	//set trigger
                    frag_trigger_set((*f_node)[cur_x][cur_y], (*f_msg)+l, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+l, frag_num, seed_len);
				/*if (abs(pre_x - cur_x) > 1)	//INS-res exist; abs(read_id) > 1? -> separate cigar
				{
					//new frag for INS-seq XXX
					_m_len = 0;
					_m_len = frag_mini_dp_multi_line(*f_node, a_msg, seed_len, (line_node){pre_x, pre_y}, (line_node){cur_x, cur_y}, _line, _line_end, 0, 0, n_seed);
					if (_m_len > 0)
					{
						for (j = 0; j < _m_len; ++j)
						{
							if (cur_num == (*line_m))
							{
								(*f_msg) = (frag_msg*)realloc((*f_msg), ((*line_m)+1) * sizeof(frag_msg));
								frag_init_msg((*f_msg)+(*line_m), (*f_msg)->frag_max); 
								if ((*f_msg) == NULL) { fprintf(stderr, "[frag dp path] Not enough memory.(line_m: %d)\n", (*line_m)+1); exit(0); }
								(*line_m)++;
							}
							frag_copy_msg((*f_msg), (*f_msg)+cur_num);
							frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, *f_node, _line+_line_end[j], _line_end[j+1]-_line_end[j]);
							((*f_msg)+cur_num)->frag_left_bound = a_msg[pre_x].read_id;
							((*f_msg)+cur_num)->frag_right_bound = a_msg[cur_x].read_id;
							cur_num++;
						}
					}
				}*/
			}
			//DEL
			else if ((*f_node)[cur_x][cur_y].match_flag == F_DELETE)
			{
				frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+l, frag_num, seed_len);
				if ((*f_node)[cur_x][cur_y].next_trigger_n > 0 || (*f_node)[cur_x][cur_y].pre_trigger_n > 0)	//set trigger
                    frag_trigger_set((*f_node)[cur_x][cur_y], (*f_msg)+l, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+l, frag_num, seed_len);
			}
			//MIS	XXX
			else if ((*f_node)[cur_x][cur_y].match_flag == F_MISMATCH || (*f_node)[cur_x][cur_y].match_flag == F_LONG_MISMATCH)
			{
				//XXX take mis-match as NEW-FLAG case
				//frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
				frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+l, frag_num, seed_len);
				if ((*f_node)[cur_x][cur_y].next_trigger_n > 0 || (*f_node)[cur_x][cur_y].pre_trigger_n > 0)	//set trigger
                    frag_trigger_set((*f_node)[cur_x][cur_y], (*f_msg)+l, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+l, frag_num, seed_len);
				//XXX find new flag in mis-match seeds
				//dp-line for INV/TRS
				/*
                    _m_len = 0;
                    _m_len = frag_mini_dp_line((*f_node), a_msg, seed_len, (line_node){pre_x, pre_y}, (line_node){cur_x, cur_y}, _line, 0, 0, n_seed);
                    if (_m_len != 0)
                    {
                        //new frag for INV/TRS
                        if (cur_num == (*line_m))
                        {
                            (*f_msg) = (frag_msg*)realloc((*f_msg), ((*line_m)+1) * sizeof(frag_msg));
                            frag_init_msg((*f_msg)+(*line_m), (*f_msg)->frag_max); 
                            if ((*f_msg) == NULL) { fprintf(stderr, "[frag dp path] Not enough memory.\n"); exit(0); }
                            (*line_m)++;
                        }
                        frag_copy_msg((*f_msg), (*f_msg)+cur_num);
                        frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, *f_node, _line, _m_len);
                        ((*f_msg)+cur_num)->frag_left_bound = a_msg[pre_x].read_id;	//bound XXX
                        ((*f_msg)+cur_num)->frag_right_bound = a_msg[cur_x].read_id;
                        cur_num++;
                    }
                */
			}
			else if ((*f_node)[cur_x][cur_y].match_flag == F_MATCH)	//: MATCH
			{
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+l, frag_num, seed_len);
			}
			else {fprintf(stderr, "[frag dp path] Error: Unknown flag, \"%d\"", (*f_node)[cur_x][cur_y].match_flag); exit(0);}
		}
		//last start
		cur_x = line[l][0].x; cur_y = line[l][0].y;
		frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+l, frag_num, seed_len);
        if ((*f_node)[cur_x][cur_y].next_trigger_n > 0 || (*f_node)[cur_x][cur_y].pre_trigger_n > 0)	//set trigger
            frag_trigger_set((*f_node)[cur_x][cur_y], (*f_msg)+l, frag_num);
        left_bound = ((line[l][line_end[l]].x == -1) ? 0 : (a_msg[line[l][line_end[l]].x].read_id));
        //left_bound = 0;
		//MIS-MATCH
		/*if (cur_x > 0)
		{
			_m_len = 0;
			_m_len = frag_mini_dp_line((*f_node), a_msg, seed_len, (line_node){-1, 0}, (line_node){cur_x, cur_y}, _line, 0, 0, n_seed);
			if (_m_len != 0)
			{
				if (cur_num == (*line_m))
				{
					(*f_msg) = (frag_msg*)realloc(*f_msg, ((*line_m)+1) * sizeof(frag_msg));
					frag_init_msg((*f_msg)+(*line_m), (*f_msg)->frag_max); 
					if ((*f_msg) == NULL) { fprintf(stderr, "[frag_dp_path] Not enough memory.(line_m: %d)\n", (*line_m)+1); exit(0); }
					(*line_m)++;
				}
				frag_copy_msg(*f_msg, (*f_msg)+cur_num);
				frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, *f_node, _line, _m_len);
				((*f_msg)+cur_num)->frag_left_bound = 0;
				((*f_msg)+cur_num)->frag_right_bound = a_msg[cur_x].read_id;
				cur_num++;
				left_bound = a_msg[_line[_m_len-1].x].read_id;
			}
		}*/
		((*f_msg)+l)->frag_left_bound = left_bound;
	}
	
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

#ifdef PLAIN_IN
int frag_cluster(const char *read_prefix, char *seed_result, seed_msg *s_msg, int seed_len, int aln_type, bntseq_t *bns, uint8_t *pac)
{
    FILE *result_p; char readline[1024];
    int n_read/*start from 1*/, n_seed, i, j; char srand;
    int read_id, chr, edit_dis; long long offset; char cigar[1024];

    sam_msg *m_msg;
    samfile_t *samf = 0;

    aln_msg *a_msg; 
    frag_msg *f_msg; int line_n, *line_tri, line_m;

    gzFile readfp; kseq_t *read_seq_t; char *read_seq;

    //alloc mem and initialization
    m_msg = sam_init_msg();
    a_msg = aln_init_msg(s_msg->seed_max);

    line_node **line = (line_node**)malloc(s_msg->seed_max * sizeof(line_node*));
    int *line_end = (int*)malloc((s_msg->seed_max) * sizeof(int));
    line_tri = (int*)malloc((s_msg->seed_max) * sizeof(int));
    line_node **_line = (line_node**)malloc(s_msg->seed_max * sizeof(line_node*));
    int *_line_end = (int*)malloc((s_msg->seed_max) * sizeof(int));
    frag_dp_node ***f_node = (frag_dp_node***)malloc(sizeof(frag_dp_node**));
    (*f_node) = (frag_dp_node**)malloc(s_msg->seed_max * sizeof(frag_dp_node*));
    for (i = 0; i < s_msg->seed_max; ++i)
    {
        line[i] = (line_node*)malloc((s_msg->seed_max+3) * sizeof(line_node));
        _line[i] = (line_node*)malloc((s_msg->seed_max+3) * sizeof(line_node));
        (*f_node)[i] = (frag_dp_node*)malloc(PER_ALN_N * sizeof(frag_dp_node)); 
        for (j = 0; j < PER_ALN_N; ++j)
        {
            (*f_node)[i][j].trigger_m = 100;
            (*f_node)[i][j].trigger = (int*)malloc(100 * sizeof(int));
        }
    }
    if (line == NULL || _line == NULL || f_node == NULL) { fprintf(stderr, "[frag_dp_path] Not enougy memory.\n"); exit(0); }

    //XXX 
    f_msg = (frag_msg*)malloc(sizeof(frag_msg));	
    frag_init_msg(&f_msg[0], s_msg->seed_max);
    line_m = 1;	//size
    line_n = 0;	//line num

    readfp = gzopen(read_prefix, "r");
    read_seq_t = kseq_init(readfp);

    //alloc mem for hash mapping
    uint32_t *hash_num;
    uint64_t **hash_node;
    int key_len = 2;
    int hash_size = (int)pow(NT_N, key_len);
    hash_num = (uint32_t*)malloc(hash_size * sizeof(uint32_t));	//16 = pow(4, 2)
    hash_node = (uint64_t**)malloc(hash_size * sizeof(uint64_t*));

    if ((result_p = fopen(seed_result, "r")) == NULL) {
        fprintf(stderr, "[lsat_aln] Can't open seed result file %s.\n", seed_result); 
        exit(-1); 
    }

    n_read = 0;
    n_seed = 0;
    int multi_aln = 1, last_id = 0, REPEAT = 0, FLAG=0;

    //while (sam_parse(
    //get seed msg of every read
    while (fgets(readline, 1024, result_p) != NULL) {
        //sam_parse();
        sscanf(readline, "%d %d %lld %c %d %s", &read_id, &chr, &offset, &srand, &edit_dis, cigar);
        if (read_id == last_id) {		// seeds from same read
            if (++multi_aln > PER_ALN_N) {
                if (!REPEAT) {
                    n_seed--;
                    REPEAT = 1;
                }
                continue;
            } else setAmsg(a_msg, n_seed, multi_aln, read_id - s_msg->n_seed[n_read-1], chr, (int64_t)offset, srand, cigar); } else {		//get a new seed 
                REPEAT = 0;
                if (read_id > s_msg->n_seed[n_read]) {	//new read
                    if (last_id != 0) {
                        //if (find_path(a_msg, n_seed, line, path_end, path, price_n, seed_len, bns->n_seqs))
                        //if (frag_find_path(a_msg, n_seed, seed_len, line, line_end, f_node, f_msg))
                        f_msg[0].last_len = s_msg->last_len[n_read]; f_msg[0].seed_all = s_msg->n_seed[n_read]-s_msg->n_seed[n_read-1];
                        if (frag_dp_path(a_msg, n_seed, seed_len, &f_msg, &line_n, line_tri, &line_m, line, line_end, f_node, _line, _line_end)) 
                            /* SW-extenging */
                            frag_check(s_msg->read_name[n_read], bns, pac, read_prefix, read_seq, s_msg->read_len[n_read], &f_msg, line_n, line_tri, a_msg, &hash_num, &hash_node, seed_len);
                    }
                    n_seed = 0;
                    while (s_msg->n_seed[n_read] < read_id) {
                        if (FLAG == 0) FLAG = 1;
                        ++n_read;
                        if (kseq_read(read_seq_t) < 0) {
                            fprintf(stderr, "[lsat_aln] Read file ERROR.\n");
                            exit(-1);
                        }
                        read_seq = read_seq_t->seq.s;
                    }
                    FLAG = 0;
                }
                multi_aln = 1;
                last_id = read_id;
                if (n_seed >= s_msg->seed_max)	{
                    fprintf(stderr, "[lsat_lan] bug: n_seed > seed_max\n");
                    exit(-1);
                }
                ++n_seed;
                setAmsg(a_msg, n_seed, multi_aln, read_id-s_msg->n_seed[n_read-1], chr, offset, srand, cigar);
            }
    }
    //if (find_path(a_msg, n_seed, line, path_end, path, price_n, seed_len, bns->n_seqs))
    //if (frag_find_path(a_msg, n_seed, seed_len, line, line_end, f_node, f_msg))
    f_msg[0].seed_all = s_msg->n_seed[n_read] - s_msg->n_seed[n_read-1];
    f_msg[0].last_len = s_msg->last_len[n_read];
    if (frag_dp_path(a_msg, n_seed, seed_len, &f_msg, &line_n, line_tri,  &line_m, line, line_end, f_node, _line, _line_end))
        /* SW-extenging */
        frag_check(s_msg->read_name[n_read], bns, pac, read_prefix, read_seq, s_msg->read_len[n_read], &f_msg, line_n, line_tri, a_msg, &hash_num, &hash_node, seed_len);

    //free variables and close file handles
    fclose(result_p);
    for (i = 0; i < m_msg->sam_m; ++i) {
        if (m_msg->sam[i].cigar_s->s) free(m_msg->sam[i].cigar_s->s);
        free(m_msg->sam[i].cigar_s);
    } 
    free(m_msg->sam); free(m_msg);
    aln_free_msg(a_msg, s_msg->seed_max);
    for (i = 0; i < s_msg->seed_max; ++i) {
        for (j = 0; j < PER_ALN_N; ++j)
            free((*f_node)[i][j].trigger);
        free((*f_node)[i]); free(line[i]); free(_line[i]);
    }
    free(*f_node); free(f_node); free(line); free(_line); free(line_end); free(line_tri); free(_line_end);
    //hash map
    free(hash_num); 
    for (i = 0; i < hash_size; ++i) free(hash_node[i]);
    free(hash_node);

    frag_free_msg(f_msg, line_m);
    gzclose(readfp); kseq_destroy(read_seq_t);
    samclose(samf);

    return 0;
}
#endif

#ifdef SAM_IN
int frag_cluster(const char *read_prefix, char *seed_result, seed_msg *s_msg, int seed_len, int aln_type, bntseq_t *bns, uint8_t *pac)
{
    int n_read/*start from 1*/, n_seed, i, j;
    int read_id;

    sam_msg *m_msg;
    samfile_t *samf = 0;

    aln_msg *a_msg; 
    frag_msg *f_msg; int line_n, *line_tri, line_m;

    gzFile readfp; kseq_t *read_seq_t; char *read_seq;

    //alloc mem and initialization
    m_msg = sam_init_msg();
    a_msg = aln_init_msg(s_msg->seed_max);

    line_node **line = (line_node**)malloc(s_msg->seed_max * sizeof(line_node*));
    int *line_end = (int*)malloc((s_msg->seed_max) * sizeof(int));
    line_tri = (int*)calloc((s_msg->seed_max), sizeof(int));
    line_node **_line = (line_node**)malloc(s_msg->seed_max * sizeof(line_node*));
    int *_line_end = (int*)malloc((s_msg->seed_max) * sizeof(int));
    frag_dp_node ***f_node = (frag_dp_node***)malloc(sizeof(frag_dp_node**));
    (*f_node) = (frag_dp_node**)malloc((s_msg->seed_max + 2) * sizeof(frag_dp_node*));
    for (i = 0; i < s_msg->seed_max; ++i)
    {
        line[i] = (line_node*)malloc((s_msg->seed_max+3) * sizeof(line_node));
        _line[i] = (line_node*)malloc((s_msg->seed_max+3) * sizeof(line_node));
    }

    for (i = 0; i < s_msg->seed_max+2; ++i)
    {
        (*f_node)[i] = (frag_dp_node*)malloc(PER_ALN_N * sizeof(frag_dp_node)); 
        for (j = 0; j < PER_ALN_N; ++j)
        {
            (*f_node)[i][j].trigger_m = 100;
            (*f_node)[i][j].trigger = (int*)calloc(100, sizeof(int));
            (*f_node)[i][j].next_trigger_n = 0;
            (*f_node)[i][j].pre_trigger_n = 0;
        }
    }
    if (line == NULL || _line == NULL || f_node == NULL) { fprintf(stderr, "[frag_dp_path] Not enougy memory.\n"); exit(0); }

    //XXX 
    f_msg = (frag_msg*)malloc(sizeof(frag_msg));	
    frag_init_msg(&f_msg[0], s_msg->seed_max);
    line_m = 1;	//size
    line_n = 0;	//line num

    readfp = gzopen(read_prefix, "r");
    read_seq_t = kseq_init(readfp);

    //alloc mem for hash mapping
    uint32_t *hash_num;
    uint64_t **hash_node;
    int key_len = 2;
    int hash_size = (int)pow(NT_N, key_len);
    hash_num = (uint32_t*)calloc(hash_size, sizeof(uint32_t));	//16 = pow(4, 2)
    hash_node = (uint64_t**)calloc(hash_size, sizeof(uint64_t*));

    if ((samf = samopen(seed_result, "r")) == 0) {
        fprintf(stderr, "[lsat_aln] Can't open seed result sam file %s.\n", seed_result);
        exit(-1);
    }
    if (samf->header == 0) {
        fprintf(stderr, "[lsat_aln] Can't read the header of result sam file %s.\n", seed_result);
        exit(-1);
    }
    int r;
    read_id = n_seed = 0;
    n_read = 1;

    while ((r = sam_read1(samf->x.tamr, samf->header, m_msg, PER_ALN_N)) >= 0) { // get seed msg of every read
        ++read_id;

        if (r > 0) {
            ++n_seed;
            for (i = 0; i < m_msg->sam_n; ++i) setAmsg(a_msg, n_seed, i+1, read_id - s_msg->n_seed[n_read-1], m_msg->sam[i].chr, m_msg->sam[i].offset, m_msg->sam[i].nsrand, m_msg->sam[i].cigar_s->s);
        }
        if (read_id == s_msg->n_seed[n_read]) { // get a whole-read
            if (kseq_read(read_seq_t) < 0) { fprintf(stderr, "[lsat_aln] Read file ERROR.\n"); exit(-1); }
            // XXX no seed is 'setAmsg'ed
            read_seq = read_seq_t->seq.s;
            f_msg[0].last_len = s_msg->last_len[n_read]; f_msg[0].seed_all = s_msg->n_seed[n_read]-s_msg->n_seed[n_read-1];

            //if (s_msg->read_name[n_read][3] == 'Y')
            //{
                if (frag_dp_path(a_msg, n_seed, seed_len, &f_msg, &line_n, line_tri, &line_m, line, line_end, f_node, _line, _line_end, aln_type))
                    frag_check(s_msg->read_name[n_read], bns, pac, read_prefix, read_seq, s_msg->read_len[n_read], &f_msg, line_n, line_tri, a_msg, &hash_num, &hash_node, seed_len);
                else
                    lsat_unmap(s_msg->read_name[n_read]);
            //}

            n_seed = 0;
            ++n_read;
        }
    }
    //free variables and close file handles
    for (i = 0; i < m_msg->sam_m; ++i) {
        if (m_msg->sam[i].cigar_s->s) free(m_msg->sam[i].cigar_s->s);
        free(m_msg->sam[i].cigar_s);
    } 
    free(m_msg->sam); free(m_msg);
    aln_free_msg(a_msg, s_msg->seed_max);
    for (i = 0; i < s_msg->seed_max+2; ++i) {
        for (j = 0; j < PER_ALN_N; ++j)
            free((*f_node)[i][j].trigger);
        free((*f_node)[i]); 
    }

    for (i = 0; i < s_msg->seed_max; ++i) {
        free(line[i]); free(_line[i]);
    }
    free(*f_node); free(f_node); free(line); free(_line); free(line_end); free(line_tri); free(_line_end);
    //hash map
    free(hash_num); 
    for (i = 0; i < hash_size; ++i) free(hash_node[i]);
    free(hash_node);

    frag_free_msg(f_msg, line_m);
    gzclose(readfp); kseq_destroy(read_seq_t);
    samclose(samf);

    return 0;
}
#endif

/* relative path convert for soap2-dp */
void relat_path(const char *ref_path, const char *soap_dir, char *relat_ref_path)	
{
    int i;
    char lsat_dir[1024], abs_soap_dir[1024], abs_ref_path[1024], ref_dir[1024], ref_file[1024];

    if (getcwd(lsat_dir, 1024) == NULL) { perror("getcwd error"); exit(-1); } 
    if (chdir(soap_dir) != 0) { perror("Wrong soap2-dp path"); exit(-1); }
    if (getcwd(abs_soap_dir, 1024) == NULL) { perror("getcwd error"); exit(-1); }

    //printf("ref: %s\n",ref_path);
    if (ref_path[0] == '.')
    {
        if (chdir(lsat_dir) != 0) { perror("Wrong soap2-dp path"); exit(-1); }
        strcpy(ref_dir, ref_path);
        for (i = strlen(ref_path)-1; i >= 0; i--)
            if (ref_path[i] == '/') { ref_dir[i] = '\0'; strncpy(ref_file, ref_path+i, 1024); }
        if (chdir(ref_dir) != 0) { perror("Wrong soap2-dp path"); exit(-1); }
        if (getcwd(abs_ref_path, 1024) == NULL) { perror("getcwd error"); exit(-1); }
        strcat(abs_ref_path, ref_file);
    }
    else
        strcpy(abs_ref_path, ref_path);

    if (chdir(lsat_dir) != 0) { perror("chdir error"); exit(-1); }
    //printf("soap: %s\nref: %s\n", abs_soap_dir, abs_ref_path);
    int dif=-1;
    for (i = 0; i < strlen(abs_soap_dir); ++i)
    {
        if (abs_soap_dir[i] != abs_ref_path[i]) break;
        if (abs_soap_dir[i] == '/') dif = i;
    }
    //printf("i: %d dif: %d\n", i, dif);
    if(dif == -1)
    {
        fprintf(stderr, "[lsat_aln] dir bug\n");
        exit(-1);
    }
    strcpy(relat_ref_path, "./");
    for (i = dif; i < strlen(abs_soap_dir); ++i)
    {
        if (abs_soap_dir[i] == '/')
            strcat(relat_ref_path, "../");
    }
    strcat(relat_ref_path, abs_ref_path+dif+1);
}

int lsat_soap2_dp(const char *ref_prefix, const char *read_prefix, char *opt_m)
{
    char relat_ref_path[1024], relat_read_path[1024];
    char lsat_dir[1024];

    if (getcwd(lsat_dir, 1024) == NULL) { perror("getcwd error"); exit(-1); } 
    relat_path(ref_prefix, SOAP2_DP_DIR, relat_ref_path);
    relat_path(read_prefix, SOAP2_DP_DIR, relat_read_path);
    if (chdir(SOAP2_DP_DIR) != 0) { perror("Wrong soap2-dp dir"); exit(-1); }

    char soap2_dp_cmd[1024];
    sprintf(soap2_dp_cmd, "./soap2-dp single %s.index %s.seed -h 2 %s > %s.seed.aln", relat_ref_path, relat_read_path, opt_m, relat_read_path);
    fprintf(stderr, "[lsat_aln] Executing soap2-dp ... ");
    if (system (soap2_dp_cmd) != 0)
        exit(-1);
    fprintf(stderr, "done.\n");

    if (chdir(lsat_dir) != 0) { perror("chdir error"); exit(-1); }
    return 0;
}

int lsat_aln_core(const char *ref_prefix, const char *read_prefix, int seed_info, int no_soap2_dp, char *seed_result, char *opt_m, int opt_l, int aln_type)
{
    seed_msg *s_msg;
    bntseq_t *bns;
    int seed_len;

    /* split-seeding */
    s_msg = seed_init_msg();
    if (seed_info)
        split_seed_info(read_prefix, s_msg, &seed_len);
    else
    {
        seed_len = opt_l;
        split_seed(read_prefix, s_msg, opt_l);
    }

    if (!strcmp(seed_result, ""))
    {
        strcpy(seed_result, read_prefix);
#ifdef PLAIN_IN
        strcat(seed_result, ".seed.out.0");
#endif
#ifdef SAM_IN
#ifdef BWA_SEED
        strcat(seed_result, ".seed.bwa.sam");
#else
        strcat(seed_result, ".seed.out.0");
#endif
#endif
    }
    //excute soap2-dp program
    if (!no_soap2_dp) lsat_soap2_dp(ref_prefix, read_prefix, opt_m);

    /* frag-clustering */
    /* SW-extending for per-frag */
    fprintf(stderr, "[lsat_aln] Restoring ref-indices ... ");
    bns = bns_restore(ref_prefix);
    uint8_t *pac = (uint8_t*)calloc(bns->l_pac/4+1, 1);
    fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);	fprintf(stderr, "done.\n");
    fprintf(stderr, "[lsat_aln] Clustering frag ... ");
    frag_cluster(read_prefix, seed_result, s_msg, opt_l, aln_type, bns, pac);	fprintf(stderr, "done.\n");

    seed_free_msg(s_msg);
    free(pac);
    bns_destroy(bns);
    return 0;
}

int lsat_aln(int argc, char *argv[])
{
    int c;
    char *ref, *read;
    // parameters
    int no_soap2_dp=0, seed_info=0, opt_l=0, aln_type=1;
    char seed_result_f[1024]="", opt_m[100], aln_resutl_f[1024]="";

    opt_l = SEED_LEN;
    strcpy(opt_m, "-m 3e ");

    while ((c =getopt(argc, argv, "nsa:m:l:t:o:")) >= 0)
    {
        switch (c)
        {
            case 'n':
                no_soap2_dp = 1;
                break;
            case 's':
                seed_info = 1;
                break;
            case 'a':
                strcpy(seed_result_f, optarg);	//soap2-dp alignment result
                break;
            case 'm':
                sprintf(opt_m, "-m %s ", optarg);
                break;
            case 'l':
                opt_l = atoi(optarg);
                break;
            case 't':
                aln_type = atoi(optarg);
                if (aln_type < 0 || aln_type > 2)
                    return usage();
                break;
            case 'o':
                strcpy(aln_resutl_f, optarg);
                break;
            default:
                return usage();
        }
    }
    if (argc - optind != 2)
        return usage();

    ref = strdup(argv[optind]);
    read =strdup(argv[optind+1]);

    lsat_aln_core(ref, read, seed_info, no_soap2_dp, seed_result_f, opt_m, opt_l, aln_type);

    free(ref); free(read);
    return 0;
}
