#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>
#include <ctype.h>
#include "lsat_aln.h"
#include "bntseq.h"
#include "frag_check.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

//char FLAG[4][10] = {"MATCH", "INSERT", "DELETION", "CHR_DIF"};

//int price[SEED_INIT_MAX][PER_ALN_N];    //p[i] = max{p[i-1][j] + cost[j,k]}
//path_msg path[SEED_INIT_MAX][PER_ALN_N];

int usage(void )		//aln usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   lsat aln [options] <ref_prefix> <in.fa/fq>\n\n");
	fprintf(stderr, "Options: -n              Do NOT excute soap2-dp program, when soap2-dp result existed already.\n");
	fprintf(stderr, "         -a [STR]        The soap2-dp alignment result. When '-n' is used. [Def=\"seed_prefix.out.0\"]\n");
	fprintf(stderr, "         -m [INT][STR]  The soap2-dp option. Maximun #errors allowed. [Def=3e]\n");
	fprintf(stderr, "         -l [INT]        The length of seed. [Def=100]\n");
	fprintf(stderr, "         -o [STR]        The output file (SAM format). [Def=\"prefix_out.sam\"]\n");

	fprintf(stderr, "\n");
	return 1;
}

seed_msg *seed_init_msg(void)
{
	seed_msg *msg = (seed_msg*)malloc(sizeof(seed_msg));

	msg->n_read = 0;
	msg->n_seed = (int *)calloc(READ_INIT_MAX, sizeof(int));
	msg->seed_max = 0;

	return msg;
}

void seed_free_msg(seed_msg *msg)
{
	free(msg->n_seed);
	free(msg);
}

int split_seed(const char *prefix, seed_msg *s_msg, int seed_len)
{
	gzFile infp;
	kseq_t *seq;
	char out_f[1024], seed_head[1024], seed_seq[1024];
	FILE *outfp;
	int m_read, n_seed, *new_p, i;

	if ((infp = gzopen(prefix, "r")) == NULL)
	{
		fprintf(stderr, "[lsat_aln] Can't open read file %s\n", prefix);
		exit(-1);
	}
	seq = kseq_init(infp);

	strcpy(out_f, prefix); strcat(out_f, ".seed");
	if ((outfp = fopen(out_f, "w")) == NULL)
	{
		fprintf(stderr, "[lsat_aln] Can't open seed file %s\n", out_f);
		exit(-1);
	}
	
	seed_seq[seed_len] = '\n';
	m_read = READ_INIT_MAX;
	while (kseq_read(seq) >= 0)
	{
		n_seed = ((seq->seq.l / seed_len) + 1) >> 1;
		if (n_seed > s_msg->seed_max) s_msg->seed_max = n_seed;
		if (s_msg->n_read == m_read-1)
		{
			m_read <<= 1;
			if ((new_p = (int*)realloc(s_msg->n_seed, m_read * sizeof(int))) == NULL)
			{
				free(s_msg->n_seed);
				fprintf(stderr, "[lsat_aln] Can't allocate more memory.\n");
				exit(1);
			}
			s_msg->n_seed = new_p;
		}
		s_msg->n_read++;
		s_msg->n_seed[s_msg->n_read] = s_msg->n_seed[s_msg->n_read-1] + n_seed;
		for (i = 0; i < n_seed; i++)
		{
			sprintf(seed_head, ">%s_%d:%d\n", seq->name.s, i, i*seed_len*2);
			strncpy(seed_seq, seq->seq.s+i*seed_len*2, seed_len);
			seed_seq[seed_len+1] = '\0';
			fputs(seed_head, outfp);
			fputs(seed_seq, outfp);
		}
	}
	gzclose(infp);
	fclose(outfp);
	kseq_destroy(seq);
	
	return 0;
}

aln_msg *aln_init_msg(int seed_max)
{
	aln_msg *msg;
	int i,j;
	msg = (aln_msg*)malloc(seed_max * sizeof(aln_msg));
	for (i = 0; i < seed_max; i++)		//drop away seed whose number of alignments > PER_ALN_N
	{
		msg[i].read_id = -1;
    	msg[i].n_aln = 0;
		msg[i].at = (aln_t*)malloc(PER_ALN_N * sizeof(aln_t));
		for (j = 0; j < PER_ALN_N; j++)
		{
			msg[i].at[j].cigar = (uint32_t*)malloc(7 * sizeof(uint32_t));//XXX default value for 3-ed
			msg[i].at[j].cigar_len = 0;
			msg[i].at[j].cmax = 7;
			msg[i].at[j].bmax = 0;
		}
	}
	return msg;
}

void aln_free_msg(aln_msg *a_msg, int seed_max)	//a_msg大小为read_max个
{
	int i,j;
	for (i = 0; i < seed_max; i++)
	{
		for (j = 0; j < PER_ALN_N; j++)
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
			fprintf(stderr, "[lsat_aln] Cigar ERROR.\n");
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
			default:	fprintf(stderr, "[lsat_aln] Cigar ERROR.\n"); exit(-1); break;
		}
		if (a_msg[seed_i].at[aln_i].cigar_len == a_msg[seed_i].at[aln_i].cmax)
		{
			a_msg[seed_i].at[aln_i].cmax <<= 2 ;
			a_msg[seed_i].at[aln_i].cigar = (uint32_t*)realloc(a_msg[seed_i].at[aln_i].cigar, a_msg[seed_i].at[aln_i].cmax * sizeof(uint32_t));
		}
		a_msg[seed_i].at[aln_i].cigar[a_msg[seed_i].at[aln_i].cigar_len] = CIGAR_GEN(x, op);
		//modify variable directly OR use a auxiliary-variable
		a_msg[seed_i].at[aln_i].cigar_len++;
		s = t+1;
	}
	a_msg[seed_i].at[aln_i].len_dif = (int)(bd - bi);
	a_msg[seed_i].at[aln_i].bmax = (int)(bd > bi ? bd : bi);
}

void setAmsg(aln_msg *a_msg, int32_t read_x, int aln_y, int read_id, int chr, int offset, int srand, int edit_dis, char *cigar)
{   //read_x: (除去unmap和repeat的)read序号, aln_y: read对应的比对结果序号(从1开始)
	if (aln_y > PER_ALN_N)
	{
		fprintf(stderr, "[lsat_aln] setAmsg ERROR!\n");
		exit(0);
	}
	a_msg[read_x-1].read_id = read_id;			//from 1
	a_msg[read_x-1].at[aln_y-1].chr = chr;
	a_msg[read_x-1].at[aln_y-1].offset = offset;	//1-base
	a_msg[read_x-1].at[aln_y-1].nsrand = ((srand=='+')?1:-1);
	a_msg[read_x-1].at[aln_y-1].edit_dis = edit_dis;
	a_msg[read_x-1].n_aln = aln_y;
	setCigar(a_msg, read_x-1, aln_y-1,  cigar);
}

int get_min_dis(aln_msg *a_msg, int i, int j, int k, int* flag, int pre_pre, int seed_len)    //(i,j)对应节点，来自i-1的第k个aln
{
	int pre = (pre_pre == 1)?(i-2):(i-1);
    if (pre < 0)
    {
        printf("ERROR\n");
        exit(0);
    }

    if (a_msg[i].at[j].chr != a_msg[pre].at[k].chr || a_msg[i].at[j].nsrand != a_msg[pre].at[k].nsrand)	//不同染色体
    {
        *flag = CHR_DIF;
    	return PRICE_DIF_CHR;
    }
    int exp = a_msg[pre].at[k].offset + a_msg[pre].at[k].nsrand * (a_msg[i].read_id - a_msg[pre].read_id) * 2 * seed_len;	
    int act = a_msg[i].at[j].offset;
    //dis <= 3 || >= -3: MATCH
    //dis > 3 :DELETION
    //dis < -3:INSERT

	//dis = offset_dis + per.edit_dis + per.cigar_len
    int dis = a_msg[pre].at[k].nsrand * (act-exp) + a_msg[pre].at[k].edit_dis + a_msg[pre].at[k].cigar_len;
    
    if (dis > 10) 
        *flag = DELETION;
    else if (dis < -10) 
        *flag = INSERT;
    else
        *flag = MATCH;
    dis=(dis>0?dis:(0-dis)); //Absolute value
    return adjest(dis);
}

// XXX FIXME: path between node that are NOT adjacent.
// e.g. 10000;10200;80000;10600/80020;10800/80040 ... 
// Here, 10600 and 100800 should be selected, rather than 800200 and 800400.
/*int path_dp(aln_msg *a_msg, int n_seed, path_msg **path, int **price, int seed_len)
{
    int i,j,k;
    int tmp, min;
    int flag;
    int rev_flag = 0;   //reverse

    //Initilization of first cloumn
    for (j = 0; j < a_msg[0].n_aln; j++)
    {
        price[0][j] = 0;
    }

    for (i = 1; i < n_seed; i++)    //n_seed条seed，相邻两条连接一条路径
    {
        for (j = 0; j < a_msg[i].n_aln; j++)
        {
            min = price[i-1][0] + get_min_dis(a_msg, i, j, 0, &flag, seed_len);
            path[i][j].from = 0;
            path[i][j].flag = flag;
            for (k = 1; k < a_msg[i-1].n_aln; k++)
            {
                if ((tmp = price[i-1][k] + get_min_dis(a_msg, i, j, k, &flag, seed_len)) < min)
                {
                    min = tmp;
                    path[i][j].from = k;
                    path[i][j].flag = flag;
                }
            }
            price[i][j] = min;
        }
    }

    return 0;
}*/
int path_dp(aln_msg *a_msg, int n_seed, path_msg **path, int **price, int seed_len)
{
	int i,j,k;
	int tmp, min;
	int flag;
	int rev_flag = 0;	//reverse

	// Initilization of first cloumn
	for (j = 0; j < a_msg[0].n_aln; j++)
		price[0][j] = 0;

	for (i = 1; i < n_seed; i++)	//n_seed seeds, two adjacent seeds have a path
	{
		//for all alns, sort before search OR search directly XXX
		//sort by tag 'NM'.
		for (j = 0; j < a_msg[i].n_aln; j++)
		{
			min = price[i-1][0] + get_min_dis(a_msg, i, j, 0, &flag,0, seed_len);
			path[i][j].from = 0;
			path[i][j].flag = flag;
			path[i][j].pre_pre = 0;
			for (k = 1; k < a_msg[i-1].n_aln; k++)
			{
				if ((tmp = price[i-1][k] + get_min_dis(a_msg, i, j, k, &flag, 0, seed_len)) < min)
				{
					min = tmp;
					path[i][j].from = k;
					path[i][j].flag = flag;
					path[i][j].pre_pre = 0;
				}
			}
			if (i >= 2)
			{
				for (k = 0; k < a_msg[i-2].n_aln; k++)
				{
					if ((tmp = price[i-2][k] + PRICE_SKIP + get_min_dis(a_msg, i, j, k, &flag, 1, seed_len)) < min)
					{
						min = tmp;
						path[i][j].from = k;
						path[i][j].flag = flag;
						path[i][j].pre_pre = 1;
					}
				}
			}
			price[i][j] = min;
		}
	}
	return 0;
}

int backtrack(aln_msg* a_msg, path_msg **path, int n_seed, int **price, int seed_len, frag_msg *f_msg)  //from end to start, find every fragment's postion
{
    if (n_seed == -1)
    {
        fprintf(stdout, "frag: 0\n");    
        return 0;
    }
    //确定回溯起点
    int i, pre;
    int min_pos=0, min_score=price[n_seed][0];
    for (i = 1; i < a_msg[n_seed].n_aln; i++)
    {
        if (price[n_seed][i] < min_score)
        {
            min_pos = i;
            min_score = price[n_seed][i];
        }
    }
    //回溯，确定fragment
    int frag_num=0, pos = min_pos, flag = 0;
    //first end
	frag_set_msg(a_msg, n_seed, pos, 1, f_msg, frag_num, seed_len);
	i = n_seed;
	while (i > 0)
	{
		if (path[i][pos].flag != MATCH)
		{
			//start
			frag_set_msg(a_msg, i, pos, 0, f_msg, frag_num, seed_len);
			flag = 1;
			frag_num++;
		}
		if (path[i][pos].pre_pre) 
		{
			pre = i-2;
		}
		else pre = i-1;
		if (pre < 0)
		{
			fprintf(stderr, "[lsat_aln] Error!\n");
			exit(-1);
		}
		pos = path[i][pos].from;
		if (flag == 1)
		{
			frag_set_msg(a_msg, pre, pos, 1, f_msg, frag_num, seed_len);
			flag = 0;
		}
		else
			frag_set_msg(a_msg, pre, pos, 2, f_msg, frag_num, seed_len);
		i = pre;
	}
	/*
    for (i = n_seed; i >=0; )
    {
        if (path[i][pos].flag != MATCH)
        {
            //start
            //fprintf(stdout, "%d\t%d\n", a_msg[i].at[pos].offset, a_msg[i].read_id);    
			frag_set_msg(a_msg, i, pos, 0, f_msg, frag_num, seed_len);
            flag = 1;
            frag_num++;
        }
		if (path[i][pos].pre_pre)	i = i-2;
		else i = i-1;
        pos = path[i][pos].from;
        if (flag == 1 && i > 0)  //获取一个fragment
        {
            //上一个fragment的end
            //fprintf(stdout, "%d\t%d\t%d\t%d\t", a_msg[i-1].at[pos].chr, a_msg[i-1].nsrand[pos], a_msg[i-1].offset[pos], a_msg[i-1].read_id);
			frag_set_msg(a_msg, i, pos, 1, f_msg, frag_num, seed_len);
            flag = 0;
        }
    } */
    //last start
	frag_set_msg(a_msg, 0, pos, 0, f_msg, frag_num, seed_len);
	return 0;
}
int frag_cluster(char *ref_prefix, char *read_prefix, char *seed_result, seed_msg *s_msg, int seed_len, bntseq_t *bns, int8_t *pac)
{
	FILE *result_p;
	char line[1024];
	int n_read/*start from 1*/, n_seed, i;
	int read_id, chr, offset, srand, edit_dis;
	char cigar[1024];
	
	aln_msg *a_msg;
	int **price;
	path_msg **path;
	frag_msg *f_msg;
	gzFile readfp;
	kseq_t *read_seq_t;
	char *read_seq;

	//alloc mem and initialization
	a_msg = aln_init_msg(s_msg->seed_max);
	price = (int**)malloc(s_msg->seed_max * sizeof(int*));
	path = (path_msg**)malloc(s_msg->seed_max * sizeof(path_msg*));
	f_msg = frag_init_msg(s_msg->seed_max);
	readfp = gzopen(read_prefix, "r");
	read_seq_t = kseq_init(readfp);
	read_seq = calloc(100001, sizeof(char));

	for (i = 0; i < s_msg->seed_max; i++) 
	{
		price[i] = (int*)malloc(PER_ALN_N * sizeof(int));
		path[i] = (int*)malloc(PER_ALN_N * sizeof(path_msg));
	}
	if ((result_p = fopen(seed_result, "r")) == NULL)
	{
		fprintf(stderr, "[lsat_aln] Can't open seed result file %s.\n", seed_result);
		exit(-1);
	}

	n_read = 0;
	n_seed = 0;
	int multi_aln = 1, last_id = 0, REPEAT = 0, FLAG=0;

	//get seed msg of every read
	while (fgets(line, 1024, result_p) != NULL)
	{
		//XXX for new-out.0 add 'cigar'
		sscanf(line, "%d %d %d %c %d %s", &read_id, &chr, &offset, &srand, &edit_dis, cigar);
		if (read_id == last_id) {		// seeds from same read
			if (++multi_aln > PER_ALN_N) {
				if (!REPEAT) {
					n_seed--;
					REPEAT = 1;
				}
				continue;
			} else setAmsg(a_msg, n_seed, multi_aln, read_id - s_msg->n_seed[n_read-1], chr, offset, srand, edit_dis, cigar);
		} else {		//get a new seed
			REPEAT = 0;
			if (read_id > s_msg->n_seed[n_read]) {	//new read
				if (last_id != 0) {
					fprintf(stdout, "read %d n_seed %d\n", n_read, n_seed);
					path_dp(a_msg, n_seed, path, price, seed_len);
					backtrack(a_msg, path, n_seed-1, price, seed_len, f_msg);
					/* SW-extenging */
					frag_check(bns, pac, read_prefix, read_seq, f_msg, a_msg, seed_len);
				}
				n_seed = 0;
				while (s_msg->n_seed[n_read] < read_id) {
					if (FLAG == 0) FLAG = 1;
					else fprintf(stdout, "read %d n_seed 0\nfrag: 0\n\n", n_read);
					n_read++;
					if (kseq_read(read_seq_t) < 0)
					{
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
			n_seed++;
			setAmsg(a_msg, n_seed, multi_aln, read_id-s_msg->n_seed[n_read-1], chr, offset, srand, edit_dis, cigar);
		}
	}
	fprintf(stdout, "n_read %d n_seed %d\n", n_read, n_seed);
	path_dp(a_msg, n_seed, path, price, seed_len);
	backtrack(a_msg, path, n_seed-1, price, seed_len, f_msg);
	/* SW-extenging */
	frag_check(bns, pac, read_prefix, read_seq, f_msg, a_msg, seed_len);

	fclose(result_p);
	aln_free_msg(a_msg, s_msg->seed_max);
	for (i = 0; i < s_msg->seed_max; i++)
	{ free(price[i]); free(path[i]); }
	free(price); free(path);
	frag_free_msg(f_msg);
	gzclose(readfp);
	kseq_destroy(read_seq_t);
	//free(read_seq);

	return 0;
}

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
	for (i = 0; i < strlen(abs_soap_dir); i++)
	{
		if (abs_soap_dir[i] != abs_ref_path[i]) break;
		if (abs_soap_dir[i] == '/') dif = i;
	}
	//printf("i: %d dif: %d\n", i, dif);
	if(dif == -1)
	{
		printf("dir bug\n");
		exit(-1);
	}
	strcpy(relat_ref_path, "./");
	for (i = dif; i < strlen(abs_soap_dir); i++)
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
	if (system (soap2_dp_cmd) != 0)
		exit(-1);

	if (chdir(lsat_dir) != 0) { perror("chdir error"); exit(-1); }
	return 0;
}

int lsat_aln_core(const char *ref_prefix, const char *read_prefix, int no_soap2_dp, char *seed_result, char *opt_m, int opt_l)
{
	seed_msg *s_msg;
	bntseq_t *bns;

	/* split-seeding */
	s_msg = seed_init_msg();
	split_seed(read_prefix, s_msg, opt_l);

	if (!strcmp(seed_result, ""))
	{
		strcpy(seed_result, read_prefix);
		strcat(seed_result, ".seed.out.0");
	}
	//excute soap2-dp program
	if (!no_soap2_dp) lsat_soap2_dp(ref_prefix, read_prefix, opt_m);

	/* frag-clustering */
	/* SW-extending for per-frag */
	bns = bns_restore(ref_prefix);
	int8_t *pac;

	pac = (int8_t*)calloc(bns->l_pac/4+1, 1);
	fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
	fprintf(stdout, "[lsat_aln] frag-clustering ...\n");
	frag_cluster(ref_prefix, read_prefix, seed_result, s_msg, opt_l, bns, pac);

	seed_free_msg(s_msg);
	free(pac);
	bns_destroy(bns);
	return 0;
}

int lsat_aln(int argc, char *argv[])
{
	int c;
	int no_soap2_dp=0, opt_l=0;
	char result_f[1024]="", opt_m[100];
	char *ref, *read;
	
	opt_l = SEED_LEN;
	strcpy(opt_m, "-m 3e ");

	while ((c =getopt(argc, argv, "na:m:l:")) >= 0)
	{
		switch (c)
		{
			case 'n':
				no_soap2_dp = 1;
				break;
			case 'a':
				strcpy(result_f, optarg);	//soap2-dp alignment result
				break;
			case 'm':
				sprintf(opt_m, "-m %s ", optarg);
				break;
			case 'l':
				opt_l = atoi(optarg);
				break;
			default:
				return usage();
		}
	}
    if (argc - optind != 2)
		return usage();

	ref = strdup(argv[optind]);
	read =strdup(argv[optind+1]);

	lsat_aln_core(ref, read, no_soap2_dp, result_f, opt_m, opt_l);
	
	return 0;
}
