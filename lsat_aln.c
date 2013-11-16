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
    msg->last_len = (int *)malloc(READ_INIT_MAX * sizeof(int));
	msg->seed_max = 0;
    msg->read_max_len = 0;

	return msg;
}

void seed_free_msg(seed_msg *msg)
{
	free(msg->n_seed);
    free(msg->last_len);
	free(msg);
}

int split_seed(const char *prefix, seed_msg *s_msg, int seed_len)
{
	gzFile infp;
	kseq_t *seq;
	char out_f[1024], seed_head[1024], seed_seq[1024];
	FILE *outfp;
	int m_read, n_seed, *new_p, i;

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
	
	fprintf(stderr, "[lsat_aln] Spliting seed ... ");
	seed_seq[seed_len] = '\n';
	m_read = READ_INIT_MAX;
	while (kseq_read(seq) >= 0)
	{
		n_seed = ((seq->seq.l / seed_len) + 1) >> 1;
		if (n_seed > s_msg->seed_max) s_msg->seed_max = n_seed;
        if (seq->seq.l > s_msg->read_max_len) s_msg->read_max_len = seq->seq.l;
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
		}
		s_msg->n_read++;
		s_msg->n_seed[s_msg->n_read] = s_msg->n_seed[s_msg->n_read-1] + n_seed;
        s_msg->last_len[s_msg->n_read] = seq->seq.l - (n_seed * 2 - 1) * seed_len;
		for (i = 0; i < n_seed; i++)
		{
			sprintf(seed_head, ">%s_%d:%d\n", seq->name.s, i, i*seed_len*2);
			strncpy(seed_seq, seq->seq.s+i*seed_len*2, seed_len);
			seed_seq[seed_len+1] = '\0';
			fputs(seed_head, outfp);
			fputs(seed_seq, outfp);
		}
	}
	fprintf(stderr, "done.\n");
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
		msg[i].skip = 0;
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
		a_msg[seed_i].at[aln_i].cigar_len++;
		s = t+1;
	}
	a_msg[seed_i].at[aln_i].len_dif = (int)(bd - bi);
	a_msg[seed_i].at[aln_i].bmax = (int)(bd > bi ? bd : bi);
}

void setAmsg(aln_msg *a_msg, int32_t read_x, int aln_y, int read_id, int chr, int64_t offset, char srand, int edit_dis, char *cigar)
{   //read_x: (除去unmap和repeat的)read序号, aln_y: read对应的比对结果序号(从1开始)
	if (aln_y > PER_ALN_N) {
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

int get_dis(aln_msg *a_msg, int pre, int pre_a, int i, int j, int *flag, int seed_len)    //(i,j)对应节点，来自pre的第pre_a个aln
{
    if (pre < 0) {fprintf(stderr, "[lsat_aln] a_msg error.\n"); exit(0);}

    if (a_msg[i].at[j].chr != a_msg[pre].at[pre_a].chr || a_msg[i].at[j].nsrand != a_msg[pre].at[pre_a].nsrand)	//different chr or different srnad
    {
        *flag = CHR_DIF;
    	return PRICE_DIF_CHR;
    }
    int64_t exp = a_msg[pre].at[pre_a].offset + a_msg[pre].at[pre_a].nsrand * (a_msg[i].read_id - a_msg[pre].read_id) * 2 * seed_len;	
    int64_t act = a_msg[i].at[j].offset;
	int64_t dis = a_msg[pre].at[pre_a].nsrand * (act-exp) - ((a_msg[pre].at[pre_a].nsrand==1)?a_msg[pre].at[pre_a].len_dif:a_msg[i].at[j].len_dif);

    if (dis > 10 && dis < DEL_THD) *flag = DELETION;
    else if (dis < -10 && dis >= (0-((a_msg[i].read_id-a_msg[pre].read_id)*2-1)*seed_len)) *flag = INSERT;
    else if (dis <= 10 && dis >= -10) *flag = MATCH; 
	else *flag = UNCONNECT;

    dis=(dis>0?dis:(0-dis)); //Absolute value
	//dis = offset_dis + pre.edit_dis + pre.cigar_len
	dis += (a_msg[pre].at[pre_a].edit_dis + a_msg[pre].at[pre_a].cigar_len);
    return adjest(dis);
}

//return: 1: yes, 0: no.
int check_in(int seed_i, int line_i, int **line, int *path_end, aln_msg *a_msg, int seed_len)
{
	int flag;

	if (path_end[line_i] == 0) return 0;
	get_dis(a_msg, line[line_i][path_end[line_i]-1], 0, seed_i, 0, &flag, seed_len);
	if (flag != UNCONNECT || flag != CHR_DIF)
	{
		line[line_i][path_end[line_i]] = seed_i;
		path_end[line_i]++;
		return 1;
	}
	else return 0;
}

void new_line(int seed_i, int **line, int *path_end, int *path_n)
{
	line[*path_n][0] = seed_i;
	path_end[*path_n] = 1;
	(*path_n)++;
}

int main_line_deter(aln_msg *a_msg, int n_seed, int seed_len, int *main_line)
{
	int i, j, last_i, path_n, main_i, flag, m_len;
	int **line, *path_end;

	path_end = (int*)calloc(n_seed, sizeof(int));
	line = (int**)malloc(n_seed * sizeof(int*));
	for (i = 0; i < n_seed; i++) line[i] = (int*)malloc(n_seed * sizeof(int));
	last_i = 0;	path_n = 0; main_i = 0;

	for (i = 0; i < n_seed; i++)
	{
		if (a_msg[i].n_aln == 1)
		{
			if (check_in(i, last_i, line, path_end, a_msg, seed_len)) 
			{
				flag = 1;
				if (path_end[last_i] > path_end[main_i]) main_i = last_i;
				continue;
			}
			for (j = 0; j < path_n; j++)
			{
				if (j == last_i) continue;
				if (check_in(i, j, line, path_end, a_msg, seed_len))
				{
					flag = 1;
					last_i = j;
					if (path_end[last_i] > path_end[main_i]) main_i = last_i;
					break;
				}
			}
			if (flag==0) new_line(i, line, path_end, &path_n);
		}
	}
	main_line = line[main_i];
	for (i = 0; i < n_seed; i++)
	{
		if (i != main_i) free(line[i]);
	}
	m_len = path_end[main_i];
	free(path_end);
	return m_len; 
}

int add_path(aln_msg *a_msg, path_msg **path, int **price, int *price_n, int start, int end, int rev, int seed_len)
{
	if (start == end) return 0;
	else if (start > end) {fprintf(stderr, "[add_path] error: start > end.\n"); exit(-1);}

	int i, j, k, last_i, from, to;
	int con_flag, flag, min, dis;

	if (rev == 1) {//add seeds onto 'end' from end+1 to start 
		from = start-1;
		to = end;
	}
	else {	//rev == -1, add seeds onto 'start' from start+1 to end 
		from = end+1;
		to = start;
	}
	price[to][0] = 0; price_n[to] = 1; path[to][0].flag = PATH_END; last_i = to;
	for (i = to-rev; i != from; i=i-rev)
	{
		price_n[i] = 0;
		for (j = 0; j < a_msg[i].n_aln; j++)
		{
			flag = 0; min = -1;
			for (k = 0; k < price_n[last_i]; k++) {
				dis = price[last_i][k] + get_dis(a_msg, i, j, last_i, k, &con_flag, seed_len);
				if (con_flag != UNCONNECT && con_flag != CHR_DIF && ((min == -1)||(dis < min)))
				{
					min = dis;
					path[i][j].from.x = last_i;
					path[i][j].from.y = k;
					path[i][j].flag = con_flag;
					flag = 1;
				}
			}
			if (flag) { 
				price[i][price_n[i]] = min;
				price_n[i]++;
			}
		}
		if (price_n[i] > 0) last_i = i;
	}
	if (rev == 1)	//repair the orientation and construce a minmum path, like reverse a single linked list.
	{
		int curr_p_x, curr_p_y, curr_f_x, curr_f_y, curr_flag, pre_p_x, pre_p_y, pre_f_x, pre_f_y, pre_flag, min_i;
		for (i = start; i < to; i++)
		{
			if (price_n[i] > 0)
			{
				min = -1;
				for (j = 0; j < price_n[i]; j++)
				{
					if (min == -1 || price[i][j] < min)
					{
						min = price[i][j];
						min_i = j;
					}
				}
				
				pre_p_x = i; pre_p_y = min_i; pre_f_x = path[i][min_i].from.x; pre_f_y = path[i][min_i].from.y; pre_flag = path[i][min_i].flag;
				path[i][min_i].flag = PATH_END;
				while (pre_flag != PATH_END) {
					curr_p_x = pre_f_x; curr_p_y = pre_f_y;
					curr_f_x = path[curr_p_x][curr_p_y].from.x; curr_f_y = path[curr_p_x][curr_p_y].from.y;
					curr_flag = path[curr_p_x][curr_p_y].flag;

					path[curr_p_x][curr_p_y].from.x = pre_p_x;
					path[curr_p_x][curr_p_y].from.y = pre_p_y;
					path[curr_p_x][curr_p_y].flag = pre_flag;

					pre_p_x = curr_p_x; pre_p_y = curr_p_y;
					pre_f_x = curr_f_x; pre_f_y = curr_f_y;
					pre_flag = curr_flag;
				}
				break;
			}
		}
	}
	return 1;
}
/********************************************/
/*	1. Use uniquely aligned seeds to:		*
 *		Determine the MAIN LINE.			*
 *	2. Add other seeds to the MAIN LINE.	*/
/********************************************/      
int path_dp(aln_msg *a_msg, int n_seed, path_msg **path, int **price, int *price_n, int seed_len, int n_seq)
{
	int i, m_len;
	int *main_line=NULL;	//uniq[n_seed], uniq_len;

	//1. Determine the main line: 
	m_len = main_line_deter(a_msg, n_seed, seed_len, main_line);
	if (m_len == 0)
	{
		fprintf(stderr, "[path_dp] no seed is uniquely aligned to the ref.\n");
		exit(-1);
	}

	//2. Add other seeds to the main line.
	add_path(a_msg, path, price, price_n, 0, main_line[0], 1, seed_len);
	for (i = 0; i < m_len-1; i++)
	{
		add_path(a_msg, path, price, price_n, main_line[i], main_line[i+1], -1, seed_len);
	}
	add_path(a_msg, path, price, price_n, main_line[i], n_seed-1, -1, seed_len);

	free(main_line);
	return 0;
}

// path between node that are NOT adjacent.
// e.g. 10000;10200;80000;10600/80020;10800/80040 ... 
// Here, 10600 and 10800 should be selected, rather than 800200 and 800400.
/*
int path_dp(aln_msg *a_msg, int n_seed, path_msg **path, int **price, int seed_len, int n_seqs)
{
	int i,j,k;
	int tmp, min;
	int aln_flag;
	//if all the alignments of current seed can NOT connect with previous seed, long-skip is allowed until a connection appears.
	//NOTE: the alignments of skiped seeds should be retained.

	// Initilization of first cloumn
	for (j = 0; j < a_msg[0].n_aln; j++)
		price[0][j] = 0;
	// DP : get one or more main path (XXX maybe two is more meaningful)
	for (i = 1; i < n_seed; i++) {	//n_seed seeds, two adjacent seeds have a path
		//for all alignments, sort before search OR search directly XXX
		//sort by tag 'NM'.
		for (j = 0; j < a_msg[i].n_aln; j++) {
			//XXX modify the process of DP
			min = price[i-1][0] + get_dis(a_msg, i, j, 0, &aln_flag,0, seed_len);
			path[i][j].from = 0;
			path[i][j].flag = aln_flag;
			path[i][j].pre_pre = 0;
			for (k = 1; k < a_msg[i-1].n_aln; k++) {
				if ((tmp = price[i-1][k] + get_dis(a_msg, i, j, k, &aln_flag, 0, seed_len)) < min) {
					min = tmp;
					path[i][j].from = k;
					path[i][j].flag = aln_flag;
					path[i][j].pre_pre = 0;
				}
			}
			if (i >= 2)	{
				for (k = 0; k < a_msg[i-2].n_aln; k++) {
					if ((tmp = price[i-2][k] + PRICE_SKIP + get_dis(a_msg, i, j, k, &aln_flag, 1, seed_len)) < min)	{
						min = tmp;
						path[i][j].from = k;
						path[i][j].flag = aln_flag;
						path[i][j].pre_pre = 1;
					}
				}
			}
			price[i][j] = min;
		}
	}
	return 0;
}
*/

int backtrack(aln_msg* a_msg, path_msg **path, int n_seed, int **price, int *price_n, int seed_len, frag_msg *f_msg)  //from end to start, find every fragment's postion
{
    if (n_seed == -1)
    {
        return 1;
    }
    //Determin the start point of backtrack.
    int i, j, min, min_i;//(path[i][min_i])
    for (i = n_seed-1; i >= 0; i--)
    {
    	if (price_n[i] > 0) {
    		min = -1;
			for (j = 0; j < price_n[i]; j++)
			{
				if (min == -1 || price[i][j] < min)
				{
					min = price[i][j];
					min_i = j;
				}
			}
			break;
    	}
    }
    //backtrack from (i, min_i)
    int last_x = i, last_y = min_i, frag_num = 0, tmp, flag = 0;
    //first end
    frag_set_msg(a_msg, last_x, last_y, 1, f_msg, frag_num, seed_len);
    while (path[last_x][last_y].flag != PATH_END) {
    	if (path[last_x][last_y].flag != MATCH) {
    		//start
    		frag_set_msg(a_msg, last_x, last_y, 0, f_msg, frag_num, seed_len);
    		flag = 1;
    		frag_num++;
    	}
    	tmp = last_x;
    	last_x = path[last_x][last_y].from.x;
    	last_y = path[tmp][last_y].from.y;

    	if (flag == 1) {
    		//next end
    		frag_set_msg(a_msg, last_x, last_y, 1, f_msg, frag_num, seed_len);
    		flag = 0; 
    	} else frag_set_msg(a_msg, last_x, last_y, 2, f_msg, frag_num, seed_len);
    }
    frag_set_msg(a_msg, last_x, last_y, 0, f_msg, frag_num, seed_len);
	return 0;
}

int frag_cluster(const char *read_prefix, char *seed_result, seed_msg *s_msg, int seed_len, bntseq_t *bns, uint8_t *pac)
{
	FILE *result_p;
	char line[1024];
	int n_read/*start from 1*/, n_seed, i;
    char srand;
	int read_id, chr, edit_dis;
    long long offset;
	char cigar[1024];
	
	aln_msg *a_msg;
	int **price, *price_n;
	path_msg **path;
	frag_msg *f_msg;
	gzFile readfp;
	kseq_t *read_seq_t;
	char *read_seq;

	//alloc mem and initialization
	a_msg = aln_init_msg(s_msg->seed_max);
	price = (int**)malloc(s_msg->seed_max * sizeof(int*));
	price_n = (int*)malloc(s_msg->seed_max * sizeof(int));
	path = (path_msg**)malloc(s_msg->seed_max * sizeof(path_msg*));
	f_msg = frag_init_msg(s_msg->seed_max);
	readfp = gzopen(read_prefix, "r");
	read_seq_t = kseq_init(readfp);
	read_seq = (char*)calloc(s_msg->read_max_len+1, sizeof(char));

	for (i = 0; i < s_msg->seed_max; i++) {
		price[i] = (int*)malloc(PER_ALN_N * sizeof(int));
		path[i] = (path_msg*)malloc(PER_ALN_N * sizeof(path_msg));
	}
	if ((result_p = fopen(seed_result, "r")) == NULL) {
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
		sscanf(line, "%d %d %lld %c %d %s", &read_id, &chr, &offset, &srand, &edit_dis, cigar);
		if (read_id == last_id) {		// seeds from same read
			if (++multi_aln > PER_ALN_N) {
				if (!REPEAT) {
					n_seed--;
					REPEAT = 1;
				}
				continue;
			} else setAmsg(a_msg, n_seed, multi_aln, read_id - s_msg->n_seed[n_read-1], chr, (int64_t)offset, srand, edit_dis, cigar);
		} else {		//get a new seed
			REPEAT = 0;
			if (read_id > s_msg->n_seed[n_read]) {	//new read
				if (last_id != 0) {
					fprintf(stdout, "read %d start %d n_seed %d\n", n_read, s_msg->n_seed[n_read-1]+1, n_seed);
					path_dp(a_msg, n_seed, path, price, price_n, seed_len, bns->n_seqs);
					if (backtrack(a_msg, path, n_seed, price, price_n, seed_len, f_msg) != 1)
					/* SW-extenging */
					    frag_check(bns, pac, read_prefix, read_seq, f_msg, a_msg, seed_len, s_msg->last_len[n_read]);
				}
				n_seed = 0;
				while (s_msg->n_seed[n_read] < read_id) {
					if (FLAG == 0) FLAG = 1;
					//else fprintf(stdout, "read %d n_seed 0\nfrag: 0\n\n", n_read);
					n_read++;
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
			n_seed++;
			setAmsg(a_msg, n_seed, multi_aln, read_id-s_msg->n_seed[n_read-1], chr, offset, srand, edit_dis, cigar);
		}
	}
	fprintf(stdout, "read %d start %d n_seed %d\n", n_read, s_msg->n_seed[n_read-1]+1, n_seed);
	path_dp(a_msg, n_seed, path, price, price_n, seed_len, bns->n_seqs);
	if (backtrack(a_msg, path, n_seed, price, price_n, seed_len, f_msg) != 1)
	/* SW-extenging */
        frag_check(bns, pac, read_prefix, read_seq, f_msg, a_msg, seed_len, s_msg->last_len[n_read]);

	fclose(result_p);
	aln_free_msg(a_msg, s_msg->seed_max);
	for (i = 0; i < s_msg->seed_max; i++)
	{ free(price[i]); free(path[i]); }
	free(price); free(price_n); free(path);
	frag_free_msg(f_msg);
	gzclose(readfp);
	kseq_destroy(read_seq_t);
	free(read_seq);

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
		fprintf(stderr, "[lsat_aln] dir bug\n");
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
	fprintf(stderr, "[lsat_aln] Executing soap2-dp ... ");
	if (system (soap2_dp_cmd) != 0)
		exit(-1);
	fprintf(stderr, "done.\n");

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
	fprintf(stderr, "[lsat_aln] Restoring ref-indices ... ");
	bns = bns_restore(ref_prefix);
	uint8_t *pac = (uint8_t*)calloc(bns->l_pac/4+1, 1);
	fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);	fprintf(stderr, "done.\n");
	fprintf(stderr, "[lsat_aln] Clustering frag ... ");
	frag_cluster(read_prefix, seed_result, s_msg, opt_l, bns, pac);	fprintf(stderr, "done.\n");

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
