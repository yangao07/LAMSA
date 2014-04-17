#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>
#include <ctype.h>
#include "lsat_aln.h"
#include "bntseq.h"
#include "frag_check.h"
#include "split_mapping.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int usage(void )		//aln usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   lsat aln [options] <ref_prefix> <in.fa/fq>\n\n");
	fprintf(stderr, "Options: -n              Do NOT excute soap2-dp program, when soap2-dp result existed already.\n");
    fprintf(stderr, "         -s              Seed information file has already existed.\n");
	fprintf(stderr, "         -a [STR]        The soap2-dp alignment result. When '-n' is used. [Def=\"seed_prefix.out.0\"]\n");
	fprintf(stderr, "         -m [INT][STR]   The soap2-dp option. Maximun #errors allowed. [Def=3e]\n");
	fprintf(stderr, "         -l [INT]        The length of seed. [Def=100]\n");
	fprintf(stderr, "         -o [STR]        The output file (SAM format). [Def=\"prefix_out.sam\"]\n");

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
    int *new_p;

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

//XXX diff of order of pre and curr?
int get_dis(aln_msg *a_msg, int pre, int pre_a, int i, int j, int *flag, int seed_len)    //(i,j)对应节点，来自pre的第pre_a个aln
{
    //if (pre < 0) {fprintf(stderr, "[lsat_aln] a_msg error.\n"); exit(0);}
	if (pre == -1) { *flag = F_MATCH; return 0; }	//for node-skip

	if ((pre == i) && (pre_a == j)) {*flag = F_MATCH; return 0;}

    if (a_msg[i].at[j].chr != a_msg[pre].at[pre_a].chr || a_msg[i].at[j].nsrand != a_msg[pre].at[pre_a].nsrand)	//different chr or different srnad
    {
        *flag = F_CHR_DIF;
    	return PRICE_DIF_CHR;
    }
    int64_t exp = a_msg[pre].at[pre_a].offset + a_msg[pre].at[pre_a].nsrand * (a_msg[i].read_id - a_msg[pre].read_id) * 2 * seed_len;	
    int64_t act = a_msg[i].at[j].offset;
	int64_t dis = a_msg[pre].at[pre_a].nsrand * (act-exp) - ((a_msg[pre].at[pre_a].nsrand==1)?(a_msg[pre].at[pre_a].len_dif):(a_msg[i].at[j].len_dif));

    if (dis > 10 && dis < DEL_THD) *flag = F_DELETE;
    else if (dis < -10 && dis >= (0-((a_msg[i].read_id-a_msg[pre].read_id)*2-1)*seed_len)) *flag = F_INSERT;
    else if (dis <= 10 && dis >= -10) *flag = F_MATCH; 
	else *flag = F_UNCONNECT;

    dis=(dis>0?dis:(0-dis)); //Absolute value
	//dis = offset_dis + pre.edit_dis + pre.cigar_len
	dis += (a_msg[pre].at[pre_a].edit_dis + a_msg[pre].at[pre_a].cigar_len);
    return adjest(dis);
}

int get_abs_dis(aln_msg *a_msg, int pre, int pre_a, int i, int j, int *flag, int seed_len)    //(i,j)对应节点，来自pre的第pre_a个aln
{
	if (pre == -1 || i == -1) { *flag = F_MATCH; return 0; }	//for bound node
	if (pre == i)
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
    	return PRICE_DIF_CHR;
    }
    int64_t exp = a_msg[pre].at[pre_a].offset + a_msg[pre].at[pre_a].nsrand * (a_msg[i].read_id - a_msg[pre].read_id) * 2 * seed_len;	
    int64_t act = a_msg[i].at[j].offset;
	int64_t dis = a_msg[pre].at[pre_a].nsrand * ((a_msg[pre].read_id < a_msg[i].read_id)?(act-exp):(exp-act)) - (((a_msg[pre].at[pre_a].nsrand) * (a_msg[pre].read_id-a_msg[i].read_id) < 0)?(a_msg[pre].at[pre_a].len_dif):(a_msg[i].at[j].len_dif));

	if (dis <= 10 && dis >= -10)
	{
		//if (abs(pre-i)==1) *flag = F_MATCH;
		if (abs(a_msg[pre].read_id - a_msg[i].read_id) == 1) *flag = F_MATCH;
		else *flag = F_MISMATCH;
	}
	else if (dis > 10 && dis < DEL_THD) *flag = F_DELETE;
    else if (dis < -10 && dis >= (0-(abs(a_msg[i].read_id-a_msg[pre].read_id)*2-1)*seed_len)) *flag = F_INSERT;
	else *flag = F_UNCONNECT;
	//printf("getdis: readid: %d %d dis %d flag %d\n", a_msg[pre].read_id, a_msg[i].read_id, dis, *flag);
    dis=abs(dis);//(dis>0?dis:(0-dis)); //Absolute value
	//dis = offset_dis + pre.edit_dis + pre.cigar_len
	//dis += (a_msg[pre].at[pre_a].edit_dis + a_msg[pre].at[pre_a].cigar_len);
	dis += (a_msg[i].at[j].edit_dis + a_msg[i].at[j].cigar_len);
	//printf("dis %d ",dis); printcigar(a_msg[pre].at[pre_a].cigar, a_msg[pre].at[pre_a].cigar_len); printf("\n");

    return dis;
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

//method of main-line determination XXX
//determine the main line before finish the whole seeds XXX, eg: when the length of main line is bigger than a number.
int main_line_deter(aln_msg *a_msg, int n_seed, int seed_len, int *m_i, int **line, int *path_end)
{
    int i, j, k, flag, con_flag;
    int path_n = 0;//, max_i = 0;
    //XXX
    int path_i[500], tmp, copy_flag[500];


    //find all the small lines with F_MATCH-CONNECT nodes
    for (i = 0; i < n_seed; ++i)
    {
        if (a_msg[i].n_aln == 1)
        {
            flag = 0;
            for (j = path_n-1; j >= 0; --j)
            {
                get_abs_dis(a_msg, line[j][path_end[j]-1], 0, i, 0, &con_flag, seed_len);
                if (con_flag == F_MATCH)
                {
                    flag = 1;
                    line[j][path_end[j]] = i;
                    ++path_end[j];
					//if (path_end[j] > path_end[max_i])
					//	max_i = j;
                    if (j != (path_n-1))
                    {
                        for (k = j+1; k != path_n; ++k)
                            path_end[k] = 0;
                        path_n = j+1;
                    }
					break;
				}
				//XXX
				if (path_end[j] > 5) break;
            }
            if (flag == 0)
            {
                new_line(i, line, path_end, path_n);
                ++path_n;
            }
        }
    }
    //sort the length of every line
    for (i = 0; i < path_n; ++i)
    	path_i[i] = i;
    for (i = 0; i < path_n-1; ++i)
    {
    	for (j = i+1; j < path_n; ++j)
    	{
    		if(path_end[path_i[i]] < path_end[path_i[j]])
    		{
    			tmp = path_i[i];
    			path_i[i] = path_i[j];
    			path_i[j] = tmp;
    		}
    	}
    }
	//combine other small lines to "max_i" line
	//XXX small lines do NOT match with each other:
	//remain the long one.
	path_end[path_n] = 0;
	/*for (i = 0; i < max_i; ++i)
	{
		get_abs_dis(a_msg, line[i][path_end[i]-1], 0, line[max_i][0], 0, &con_flag, seed_len);
		if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
			copy_line(line, i, path_n, path_end);
	}
	copy_line(line, max_i, path_n, path_end);
	for (i = max_i+1; i < path_n; ++i)
	{
		get_abs_dis(a_msg, line[max_i][path_end[max_i]-1], 0, line[i][path_end[i]-1], 0, &con_flag, seed_len);
		if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
			copy_line(line, i, path_n, path_end);
	}*/
	copy_flag[path_i[0]] = 1;
	for (i = 1; i < path_n; ++i)
	{
		flag = 1;
		for (j = 0; j < i; ++j)
		{
			if (path_i[i] < path_i[j])
				get_abs_dis(a_msg, line[path_i[i]][path_end[path_i[i]]-1], 0, line[path_i[j]][0], 0, &con_flag, seed_len);
			else
				get_abs_dis(a_msg, line[path_i[j]][path_end[path_i[j]]-1], 0, line[path_i[i]][0], 0, &con_flag, seed_len);
			if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
			{
				flag = 0;
				break;
			}
		}
		copy_flag[path_i[i]] = flag;
	}
	for (i = 0; i < path_n; ++i)
	{
		if (path_end[i] == 1) continue;	//delete all the line contain only one seed.
		if (copy_flag[i] == 1) copy_line(line, i, path_n, path_end);
	}
	(*m_i) = path_n;
	return path_end[path_n];
}

//@function: use the main-line calculated before, 
int add_path(aln_msg *a_msg, path_msg **path, int *price_n, int start, int end, int rev, int seed_len)
{
	if (start == end) return 0;
	else if (start > end) {fprintf(stderr, "[add_path] error: start > end. %d %d\n", start, end); exit(-1);}

	int i, j, k, l, anchor, termi;
	int last_i, last_j, start_dis, last_dis, last_start_dis, tmp_start_dis, tmp_j, tmp_flag;
	int back_i, back_j;// for rev-path
	int con_flag, last_flag, flag, min, tmp;
	if (rev == 1) {//rev-path add seeds onto 'end' from end+1 to start 
		termi = start-1; anchor = end;
	} else {	//rev == -1/-2, add seeds onto 'start' from start+1 to end 
		termi = end+1; anchor = start;
	}
	price_n[anchor] = 1; back_i= anchor; back_j = 0;;
	if (rev == -1 || rev == 1)	//onepoint
	{
		last_i = anchor; last_j = 0;
		for (i = anchor-rev; i != termi; i=i-rev)
		{
			price_n[i] = 0;
			flag = 0; min = -1;
			for (j = 0; j < a_msg[i].n_aln; ++j)
			{
				start_dis = get_abs_dis(a_msg, anchor, 0, i, j, &con_flag, seed_len);
				if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF) {
                    last_dis = get_abs_dis(a_msg, last_i, last_j, i, j, &last_flag, seed_len);
                    if (last_flag == F_MATCH)
                    {
                        path[i][j].from.x = last_i;
                        path[i][j].from.y = last_j;
                        path[i][j].flag = F_MATCH;
                        last_i = i; last_j = j; last_start_dis = start_dis;
                        price_n[i] = j+1;
                        back_i = i; back_j = j;
                        flag = 0;
                        break;
                    }
                    else if ((min == -1) || (last_dis < min)) {
                        min = last_dis;
                        flag = 1;
                        tmp_j = j;
                        tmp_start_dis = start_dis;
                        tmp_flag = last_flag;
                    }
				}
			}
			if (flag == 1)
			{
				if (tmp_flag != F_UNCONNECT)
				{
					path[i][tmp_j].from.x = last_i;
					path[i][tmp_j].from.y = last_j;
					path[i][tmp_j].flag = tmp_flag;
					last_i = i; last_j = tmp_j; last_start_dis = tmp_start_dis;
					price_n[i] = tmp_j+1;
					back_i = i; back_j = tmp_j;
				}
				else
				{
					if (tmp_start_dis < last_start_dis)
					{
						k = path[last_i][last_j].from.x; l = path[last_i][last_j].from.y;
						while (1)
						{
							get_abs_dis(a_msg, k, l, i, tmp_j, &last_flag, seed_len);
							if (last_flag != F_UNCONNECT)
							{
								path[i][tmp_j].from.x = k;
								path[i][tmp_j].from.y = l;
								path[i][tmp_j].flag = last_flag;
								last_i = i; last_j = tmp_j; last_start_dis = tmp_start_dis;
								price_n[i] = tmp_j+1;
								back_i = i; back_j = tmp_j;
								break;
							}
							tmp = k;
							k = path[k][l].from.x; l = path[tmp][l].from.y;
						}
					}
				}
			}
        }
    }
	else //rev == -2 : two endpoint XXX
	{
		int end_flag;
		last_i = start; last_j = 0;
		for	(i = start+1; i != end+1; ++i)
		{
			price_n[i] = 0;
			flag = 0; min = -1;
			for (j = 0; j < a_msg[i].n_aln; ++j)
			{
				start_dis = get_abs_dis(a_msg, start, 0, i, j, &con_flag, seed_len);
				if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF) 
				{
					get_abs_dis(a_msg, i, j, end, 0, &end_flag, seed_len);
					if (end_flag != F_UNCONNECT && end_flag != F_CHR_DIF) {
                        last_dis = get_abs_dis(a_msg, last_i, last_j, i, j, &last_flag, seed_len);
                        //XXX aln-res are both F_MATCH, then compare on the cigar_len and edit-dis.
                        /*if (last_flag == F_MATCH) {
                            path[i][j].from.x = last_i;
                            path[i][j].from.y = last_j;
                            path[i][j].flag = F_MATCH;
                            last_i = i; last_j = j; last_start_dis = start_dis;
                            price_n[i] = j+1;
                            flag = 0;
                            break;
                        }
                        else */
                        if ((min == -1) || (last_dis < min)) {
                            min = last_dis;
                            flag = 1;
                            tmp_j = j;
                            tmp_start_dis = start_dis;
                            tmp_flag = last_flag;
                        }
                    }
                }
			}
			if (flag == 1)
			{
				if (tmp_flag != F_UNCONNECT)
				{
					//printf("i: %d, j: %d min: %d ",i,tmp_j, min); printcigar(a_msg[i].at[tmp_j].cigar, a_msg[i].at[tmp_j].cigar_len); printf("\n");
					path[i][tmp_j].from.x = last_i;
					path[i][tmp_j].from.y = last_j;
					path[i][tmp_j].flag = tmp_flag;
					last_i = i; last_j = tmp_j; last_start_dis = tmp_start_dis;
					price_n[i] = tmp_j+1;
				}
				else	//one of these two alns ([i,tmp_j], [last_i,last_j]) is wrong.
				{
					if (tmp_start_dis < last_start_dis)
					{
						k = path[last_i][last_j].from.x; l = path[last_i][last_j].from.y;
						while (1)
						{
							get_abs_dis(a_msg, k, l, i, tmp_j, &last_flag, seed_len);
							if (last_flag != F_UNCONNECT)
							{
								path[i][tmp_j].from.x = k;
								path[i][tmp_j].from.y = l;
								path[i][tmp_j].flag = last_flag;
								last_i = i; last_j = tmp_j; last_start_dis = tmp_start_dis;
								price_n[i] = tmp_j+1;
								break;
							}
							tmp = k;
							k = path[k][l].from.x; l = path[tmp][l].from.y;
						}
					}
				}
			}
		}
	}
	if (rev == 1)	//repair the orientation and construce a minmum path, like reverse a single linked list.
	{
		int curr_p_x, curr_p_y, curr_f_x, curr_f_y, curr_flag, pre_p_x, pre_p_y, pre_f_x, pre_f_y, pre_flag;
		if (price_n[back_i] > 0)
		{
			pre_p_x = back_i;    pre_p_y = back_j; 
			pre_f_x = path[back_i][back_j].from.x;
			pre_f_y = path[back_i][back_j].from.y;
			pre_flag = path[back_i][back_j].flag;
			path[back_i][back_j].flag = F_PATH_END;
			while (pre_p_x != anchor) {
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
		}
		else {fprintf(stderr, "[add_path] Bug: rev error.\n"); exit(-1);}
	}
	return 1;
}

int frag_add_seed(aln_msg *a_msg, int seed_i, int aln_i, int *path_n, line_node **line, int *line_end, int seed_len)
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

int frag_DP_init(frag_DP_node *f_node, line_node **line, int *line_end, int path_i, line_node from, aln_msg *a_msg, int seed_len)
{
	f_node[path_i].seed_num = line_end[path_i];
	f_node[path_i].n_line = 1;
	f_node[path_i].line_i[0] = path_i;
	f_node[path_i].from_i = -1;
	int con_flag;
	get_abs_dis(a_msg, from.x, from.y, line[path_i][0].x, line[path_i][0].y, &con_flag, seed_len);
	f_node[path_i].score = (con_flag == F_MATCH ? (line_end[path_i]):(line_end[path_i] - SV_PEN));
	return 0;
}

int frag_update_line(frag_DP_node *f_node, line_node **line, int *line_end, int path_i, aln_msg *a_msg, int seed_len)
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
			if (f_node[i].score + f_node[path_i].seed_num > max_score)
			{
				max_from = i;
				max_score = f_node[i].score + f_node[path_i].seed_num;
			}
		}
		else if (f_node[i].score + f_node[path_i].seed_num - SV_PEN > max_score)
		{
			max_score = f_node[i].score + f_node[path_i].seed_num - SV_PEN;
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

int mini_frag_main_line(aln_msg *a_msg, line_node left, line_node right, int seed_len, int n_seed, line_node **line, int *line_end, frag_DP_node *f_node)
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
/********************************************/
/*	1. Use uniquely aligned seeds to:		*
 *		Determine the MAIN LINE.			*
 *	2. Add other seeds to the MAIN LINE.	*/
/********************************************/      
/*int find_path(aln_msg *a_msg, int n_seed, int **line, int *path_end, path_msg **path, int *price_n, int seed_len, int n_seq)
//int find_path(aln_msg *a_msg, int n_seed, int *main_line, int **main_price, int **main_path, path_msg **path, int **price, int *price_n, int seed_len, int n_seq)
{
	int i, m_i, m_len;

	//1. Determine the main line: 
    //m_len = main_line_deter(a_msg, n_seed, seed_len, &m_i, line, path_end);
    //m_len = main_line_deter(a_msg, n_seed, seed_len, main_price, main_path, main_line);
	if (m_len == 0)
	{
		fprintf(stderr, "[find_path] no seed is uniquely aligned to the ref. %d %lld\n", a_msg[0].at[0].chr, (long long)a_msg[0].at[0].offset);
		return 0;
	}
   	//printf("main line: ");
	//for (i = 0; i < m_len; ++i)
	//		printf("%d ", a_msg[line[m_i][i]].read_id);
 	//printf("\n");
	//2. Add other seeds to the main line.
	
   
	if (add_path(a_msg, path, price_n, 0, line[m_i][0], 1, seed_len) == 0) //rev-path one endpoint
		path[line[m_i][0]][0].flag = F_PATH_END;
	for (i = 0; i < m_len-1; i++)
	{
		//fprintf(stderr, "debug: %d %d %d %lld %lld %lld \n", i, line[m_i][i], line[m_i][i+1], a_msg[line[m_i][i]].read_id, a_msg[line[m_i][i+1]].read_id, a_msg[line[m_i][i]].at[0].offset);
		add_path(a_msg, path, price_n, line[m_i][i], line[m_i][i+1], -2, seed_len);	//two endpoint
	}
	add_path(a_msg, path, price_n, line[m_i][m_len-1], n_seed-1, -1, seed_len);	//one endpoint
    
    //if (add_path(a_msg, path, price, price_n, 0, main_line[m_len-1], 1, seed_len) == 0)
    //    path[main_line[m_len-1]][0].flag = F_PATH_END;

    //for (i = m_len-1; i > 0; --i)
    //    add_path(a_msg, path, price, price_n, main_line[i], main_line[i-1], -2, seed_len);
    //add_path(a_msg, path, price, price_n, main_line[0], n_seed-1, -1, seed_len);

	return 1;
}*/

/*int frag_find_path(aln_msg *a_msg, int n_seed, int seed_len, line_node **line, int *line_end, frag_DP_node *f_node, frag_msg *f_msg)
{
	if (n_seed == -1) return 0;
	int m_i;
	m_i = frag_main_line(a_msg, n_seed, seed_len, line, line_end, f_node);
	if (line_end[m_i] == 0) return 0;
	
	int i, frag_num = 0, pre_x, pre_y, con_flag;
	int cur_x = line[m_i][line_end[m_i]-1].x, cur_y = line[m_i][line_end[m_i]-1].y;

	frag_set_msg(a_msg, cur_x, cur_y, FRAG_END, f_msg, frag_num, seed_len);
#ifdef __DEBUG__
	fprintf(stdout, "%d %d    end    %d ref %lld read %d\n", cur_x, cur_y, a_msg[cur_x].at[cur_y].nsrand, (long long)a_msg[cur_x].at[cur_y].offset, a_msg[cur_x].read_id);
#endif
	
	for (i = line_end[m_i]-1; i > 0; --i)
	{
		cur_x = line[m_i][i].x;
		cur_y = line[m_i][i].y;
		pre_x = line[m_i][i-1].x;
		pre_y = line[m_i][i-1].y;
		get_abs_dis(a_msg, pre_x, pre_y, cur_x, cur_y, &con_flag, seed_len);
		if (con_flag != F_MATCH && con_flag != F_MISMATCH)		// F_MISMATCH && MATCH XXX
		{
			frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, f_msg, frag_num, seed_len);
#ifdef __DEBUG__
			fprintf(stdout, "%d %d    start  %d ref %lld read %d\n", cur_x, cur_y, a_msg[cur_x].at[cur_y].nsrand, (long long)a_msg[cur_x].at[cur_y].offset, a_msg[cur_x].read_id);
#endif
			++frag_num;
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, f_msg, frag_num, seed_len);
#ifdef __DEBUG__
			fprintf(stdout, "%d %d    end    %d ref %lld read %d\n", pre_x, pre_y, a_msg[pre_x].at[pre_y].nsrand, (long long)a_msg[pre_x].at[pre_y].offset, a_msg[pre_x].read_id);
#endif
		}
		else
		{
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, f_msg, frag_num, seed_len);
			//fprintf(stdout, "%d %d    seed   %d ref %lld read %d\n", pre_x, pre_y, a_msg[pre_x].at[pre_y].nsrand, (long long)a_msg[pre_x].at[pre_y].offset, a_msg[pre_x].read_id);
		}
	}
	//i = 0 : START
	cur_x = line[m_i][0].x;
	cur_y = line[m_i][0].y;
	frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, f_msg, frag_num, seed_len);
#ifdef __DEBUG__
	fprintf(stdout, "%d %d    start   %d ref %lld read %d\n", cur_x, cur_y, a_msg[cur_x].at[cur_y].nsrand, (long long)a_msg[cur_x].at[cur_y].offset, a_msg[cur_x].read_id);
#endif
	return 1;
}*/

int frag_dp_init(frag_dp_node **f_node, aln_msg *a_msg, int seed_i, line_node from, int seed_len, int dp_flag)
{
	int i;
	if (from.x == -1)	//UNLIMITED
	{
		for (i = 0; i < a_msg[seed_i].n_aln; ++i)
		{
			f_node[seed_i][i] = (frag_dp_node){from, seed_i, i, 0, F_MATCH, dp_flag, 1};
		}
	}
	else
	{
		int con_flag;
		for (i = 0; i < a_msg[seed_i].n_aln; ++i)
		{
			get_abs_dis(a_msg, from.x , from.y, seed_i, i, &con_flag, seed_len);
			if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
				//f_node[seed_i][i] = (frag_dp_node){from, seed_i, i, 0, con_flag, dp_flag, 1};		XXX MATCH, MISMATCH: same score?
				f_node[seed_i][i] = (frag_dp_node){from, seed_i, i, (con_flag==F_MATCH||con_flag==F_MISMATCH)?1:0, con_flag, dp_flag, 1};
			else
				f_node[seed_i][i] = (frag_dp_node){(line_node){-1,0}, seed_i, i, -1, con_flag, 0-dp_flag, 0};
		}
	}
	return 0;
}

//pruning	XXX
int frag_dp_update(frag_dp_node **f_node, aln_msg *a_msg, int seed_i, int aln_i, int start, int seed_len, int dp_flag)
{
	int i, j, con_flag;
	line_node max_from;
	int max_score, max_flag;

	max_from = f_node[seed_i][aln_i].from;
	max_score = f_node[seed_i][aln_i].score;

	for (i = seed_i - 1; i >= start; --i)
	{
		for (j = 0; j < a_msg[i].n_aln; ++j)
		{
			if (f_node[i][j].dp_flag == dp_flag)
			{
				get_abs_dis(a_msg, i, j, seed_i, aln_i, &con_flag, seed_len);
				if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
					continue;
				//											XXX MATCH, MISMATCH: same score?
				if ((f_node[i][j].score + 1 - ((con_flag == F_MATCH||con_flag == F_MISMATCH)?0:SV_PEN)) > max_score)	
				{
					max_from = (line_node){i,j};
					max_score = f_node[i][j].score + 1 - (con_flag==F_MATCH||con_flag ==F_MISMATCH?0:SV_PEN);
					max_flag = con_flag;
				}
			}
		}
	}
	if (max_from.x != f_node[seed_i][aln_i].from.x || max_from.y != f_node[seed_i][aln_i].from.y)
	{
		f_node[seed_i][aln_i].from = max_from;
		f_node[seed_i][aln_i].score = max_score;
		f_node[seed_i][aln_i].match_flag = max_flag;
		f_node[seed_i][aln_i].node_n += f_node[max_from.x][max_from.y].node_n;
	}
	return 0;
}

//for multi-dp-lines
int frag_mini_dp_multi_line(frag_dp_node **f_node, aln_msg *a_msg, int seed_len, line_node left, line_node right, line_node *line, int *line_end, /*int _head,0,  int _tail, 0*/ int n_seed)
{
	//line_node head = ((_head)?left:(line_node){-1,0});
	line_node head; head.x = -1, head.y = 0;
	int i, j, mini_dp_flag = MULTI_FLAG;
	int l_i = 0;
	int max_score, max_n, node_i;
	line_node max_node, _right, _left;

	line_end[0] = 0;

	//first dp int
	for (i = left.x+1; i < right.x; ++i)
		frag_dp_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);
	while (1)
	{
		//dp update
		for (i = left.x+2; i < right.x; ++i)
		{
			for (j = 0; j < a_msg[i].n_aln; ++j)
			{
				if (f_node[i][j].dp_flag == mini_dp_flag)
					frag_dp_update(f_node, a_msg, i, j, left.x+1, seed_len, mini_dp_flag);
			}
		}
		//max node
		max_score = 0;
		for (i = right.x-1; i > left.x; --i)
		{
			for (j = 0; j < a_msg[i].n_aln; ++j)
			{
				if (f_node[i][j].dp_flag == mini_dp_flag && f_node[i][j].score > max_score)
				{
					max_score = f_node[i][j].score;
					max_node = (line_node){i,j};
					max_n = f_node[i][j].node_n;
				}
			}
		}
		if (max_score == 0) return l_i;
		//backtrack
		_right = max_node;
		node_i = max_n - 1;
		line_end[l_i+1] = line_end[l_i] + max_n;
		while (_right.x != head.x)
		{
			if (node_i < 0) { fprintf(stderr, "[frag mini dp] node_i < 0, BUG.\n"); exit(0); }
			line[line_end[l_i] + node_i] = _right;
			--node_i;
			_left = f_node[_right.x][_right.y].from;
			_right = _left;
		}
		if (node_i >= 0) { fprintf(stderr, "[frag mini dp] node_i >= 0, BUG.\n"); exit(0); }
		//dp init
		for (i = left.x+1; i < right.x; ++i)
			frag_dp_init(f_node, a_msg, i, head, seed_len, mini_dp_flag);	
		//remove nodes that are already in line
		for (i = 0; i < line_end[l_i+1]; ++i)
			f_node[line[i].x][line[i].y].dp_flag = 0-mini_dp_flag;
		l_i++;

	}
	//dp init
	/*for (i = left.x + 1; i < right.x; ++i)
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

	//find max-value node, backtrack
	//find max-value node in rest nodes, backtrack
	//...
	//find backtrack start node
		int max_score = 0, max_n = 0, l_i=0;
		line_node max_node = head;//left;
		//UNLIMITED 
		//if (_tail == 0)
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
		//backtrack
		line_node _right = max_node, _left;
		int node_i = max_n-1;
		//while (_right.x != left.x)
		while (_right.x != head.x)
		{
			if (node_i < 0) { fprintf(stderr, "[frag mini dp] node_i BUG 1.\n"); exit(0); }
			line[l_i][node_i--] = _right;
			_left = f_node[_right.x][_right.y].from;
			_right = _left;
		}
		if (node_i >= 0) { fprintf(stderr, "[frag mini dp] node_i BUG 2.\n"); exit(0); }*/

	/*for (i = left.x+1; i < right.x; ++i)
	{
		for (j = 0; j < a_msg[i].n_aln; ++j)
		{
			fprintf(stdout, "node:(%d %d)\t%d %d %d\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, a_msg[i].at[j].offset, f_node[i][j].from.x, f_node[i][j].from.y, f_node[i][j].score, "MXIDCRUSE"[f_node[i][j].match_flag], f_node[i][j].dp_flag, f_node[i][j].node_n);
		}
	}*/
	return max_n;
}

int frag_mini_dp_line(frag_dp_node **f_node, aln_msg *a_msg, int seed_len, line_node left, line_node right, line_node *line, int _head, int _tail, int n_seed)
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

int frag_dp_line(aln_msg *a_msg, int n_seed, int seed_len, line_node *line, frag_dp_node **f_node, line_node *_line)//_headi=0, _tail=0
{
	int i, j;
	int min_len = 1, min_exist=0;

	//dp init
	{
		for (i = 0; i < n_seed; ++i)
		{
			if (a_msg[i].n_aln <= min_len)
			{
				frag_dp_init(f_node, a_msg, i, (line_node){-1,0}, seed_len, MIN_FLAG);
				min_exist = 1;
			}
			else
				frag_dp_init(f_node, a_msg, i, (line_node){-1,0}, seed_len, MULTI_FLAG);
		}
	}
	//dp update and backtrack
	{
		if (min_exist)
		{
			//min extend
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
						if (f_node[i][j].dp_flag == MIN_FLAG)
							frag_dp_update(f_node, a_msg, i, j, 0/*update start pos*/, seed_len, MIN_FLAG);
					}
				}
			//find backtrack start node
				//_tail=0
				int max_score=0;
				line_node max_node = (line_node){-1,0};
				for (i = n_seed-1; i >= 0; --i)
				{
					for (j = 0; j < a_msg[i].n_aln; ++j)
					{
						if (f_node[i][j].dp_flag == MIN_FLAG && f_node[i][j].score > max_score)
						{
							max_score = f_node[i][j].score;
							max_node = (line_node){i,j};
						}
					}
				}
			//backtrack
				int node_i = 0, mini_len;
				line_node right, left;
				//line_node *_line = (line_node*)malloc(n_seed * sizeof(line_node));
				if (max_node.x < n_seed - 1)
				{
					mini_len = frag_mini_dp_line(f_node, a_msg, seed_len, max_node, (line_node){n_seed, 0}, _line, 1, 0, n_seed);
					for (i = mini_len-1; i >= 0; --i)
						line[node_i++] = _line[i];
				}
				right = max_node;
				while (right.x != -1)
				{
					line[node_i++] = right;
					left = f_node[right.x][right.y].from;

					if (left.x < right.x - 1 )//f_node[right.x][right.y].match_flag != F_MATCH)	//XXX if left.x < right.x-1, match_flag of right node couldn't be match?
					{
						//mini frag dp
							mini_len = frag_mini_dp_line(f_node, a_msg, seed_len, left, right, _line, 1, 1, n_seed);
							for (i = mini_len-1; i >= 0; --i)
								line[node_i++] = _line[i];
					}
					right = left;
				}
				//free(_line);
				//invert line
				line_node tmp;
				for (i = 0; i < node_i/2; ++i)
				{
					tmp =line[i]; line[i] = line[node_i-i-1]; line[node_i-i-1] = tmp;
				}

				//print
				/*for (i = 0; i < n_seed; ++i)
				{
					for (j = 0; j < a_msg[i].n_aln; ++j)
					{
						fprintf(stdout, "node:(%d %d)\t%d %d %d\tfrom:(%d %d)\tscore: %d\tM-flag:%c\tDP-flag:%d\tnode_n:%d\n", i, j, a_msg[i].at[j].nsrand, a_msg[i].at[j].chr, a_msg[i].at[j].offset, f_node[i][j].from.x, f_node[i][j].from.y, f_node[i][j].score, "MXIDCRUSE"[f_node[i][j].match_flag], f_node[i][j].dp_flag, f_node[i][j].node_n);
					}
				}*/

				return node_i;
		}
		else	//whole-multi dp, should not go to here for current data.
		{
			return 0;
		}
	}
}

int frag_mini_dp_path(aln_msg *a_msg, int seed_len, frag_msg *f_msg, frag_dp_node **f_node, line_node *line, int m_len)
{
	int i, j;
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
int frag_dp_path(aln_msg *a_msg, int n_seed, int seed_len, frag_msg **f_msg, int *line_n, int *line_m, line_node *line, frag_dp_node **f_node, line_node *_line)
{
	int i, j;
	//DP
	/*line_node *line = (line_node*)malloc(n_seed * sizeof(line_node));
	frag_dp_node **f_node = (frag_dp_node**)malloc((n_seed) * sizeof(frag_dp_node*));
	for (i = 0; i < n_seed; ++i)
		f_node[i] = (frag_dp_node*)malloc(a_msg[i].n_aln * sizeof(frag_dp_node)); 
	if (line == NULL || f_node == NULL) { fprintf(stderr, "[frag_dp_path] Not enougy memory.\n"); exit(0); }*/

	int m_len = frag_dp_line(a_msg, n_seed, seed_len, line, f_node, _line);
	if (m_len == 0) return 0;

	int frag_num = 0;
	int cur_x , cur_y , pre_x=line[m_len-1].x, pre_y=line[m_len-1].y;

	//first end
	/*frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg), frag_num, seed_len);
	for (i = m_len-1; i > 0; --i)
	{
		cur_x = pre_x; cur_y = pre_y;
		pre_x = line[i-1].x; pre_y = line[i-1].y;
		if (f_node[cur_x][cur_y].match_flag != F_MATCH && f_node[cur_x][cur_y].match_flag != F_MISMATCH)
		{
			frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg), frag_num, seed_len);
			++frag_num;
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg), frag_num, seed_len);
		}
		else	//seed
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg), frag_num, seed_len);
	}
	cur_x = line[0].x; cur_y = line[0].y;
	frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg), frag_num, seed_len);*/

	//for new frags
	//  cur_line: cur line index, cur_num:  whole lines, *line_num: mem of *f_msg
	int cur_line=0, cur_num=1, _m_len;
	//line_node *_line = (line_node*)malloc(n_seed * sizeof(line_node));


	//right bound
	int right_bound = (*f_msg)[0].seed_all+1;
	//MIS-MATCH
	if (abs(n_seed-pre_x) > 1)
	{
		_m_len = 0;
		_m_len = frag_mini_dp_line(f_node, a_msg, seed_len, (line_node){pre_x, pre_y}, (line_node){n_seed, 0}, _line, 0, 0, n_seed);
		if (_m_len != 0)
		{
			if (cur_num == (*line_m))
			{
				(*f_msg) = (frag_msg*)realloc(*f_msg, ((*line_m)+1) * sizeof(frag_msg));
				if ((*f_msg) == NULL) { fprintf(stderr, "[frag_dp_path] Not enough memory.(line_m: %d)\n", (*line_m)+1); exit(0); }
				(*line_m)++;
			}
			frag_copy_msg(*f_msg, (*f_msg)+cur_num);
			frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, f_node, _line, _m_len);
			((*f_msg)+cur_num)->frag_left_bound = a_msg[pre_x].read_id;
			((*f_msg)+cur_num)->frag_right_bound = (*f_msg)[0].seed_all+1;
			cur_num++;
			right_bound = a_msg[_line[0].x].read_id;
		}
	}
	//first end
	frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
	((*f_msg)+cur_line)->frag_right_bound = right_bound;//(*f_msg)[0].seed_all+1;	//read_id
	for (i = m_len-1; i > 0; --i)
	{
		cur_x = pre_x; cur_y = pre_y;
		pre_x = line[i-1].x; pre_y = line[i-1].y;
		//INS
		if (f_node[cur_x][cur_y].match_flag == F_INSERT)
		{
			frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
			++frag_num;
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
			if (abs(pre_x - cur_x) > 1)	//INS-res exist; abs(read_id) > 1? -> separate cigar
			{
				//new frag for INS-seq XXX
				_m_len = 0;
				int *_line_end = (int*)malloc(n_seed * sizeof(line_node));
				_m_len = frag_mini_dp_multi_line(f_node, a_msg, seed_len, (line_node){pre_x, pre_y}, (line_node){cur_x, cur_y}, _line, _line_end, n_seed);
				if (_m_len > 0)
				{
					for (j = 0; j < _m_len; ++j)
					{
						if (cur_num == (*line_m))
						{
							(*f_msg) = (frag_msg*)realloc((*f_msg), ((*line_m)+1) * sizeof(frag_msg));
							if ((*f_msg) == NULL) { fprintf(stderr, "[frag dp path] Not enough memory.(line_m: %d)\n", (*line_m)+1); exit(0); }
							(*line_m)++;
						}
						frag_copy_msg((*f_msg), (*f_msg)+cur_num);
						frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, f_node, _line+_line_end[j], _line_end[j+1]-_line_end[j]);
						((*f_msg)+cur_num)->frag_left_bound = a_msg[pre_x].read_id;
						((*f_msg)+cur_num)->frag_right_bound = a_msg[cur_x].read_id;
						cur_num++;
					}
				}
				free(_line_end);
			}
		}
		//DEL
		else if (f_node[cur_x][cur_y].match_flag == F_DELETE)
		{
			frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
			++frag_num;
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
		}
		//MIS	XXX
		/*else if (f_node[cur_x][cur_y].match_flag == F_MISMATCH)//: mis-match/ref-N/INV/TRS
		{
			frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
			++frag_num;
			//check for INV/TRS
			if (abs(pre_x - cur_x) > 1)
			{
				//dp-line
				_m_len = 0;
				_m_len = frag_mini_dp_line(f_node, a_msg, seed_len, (line_node){pre_x, pre_y}, (line_node){cur_x, cur_y}, _line, 0, 0, n_seed);
				//new frag for INV/TRS-seq
				if (_m_len != 0)
				{
					if (cur_num == (*line_m))
					{
						(*f_msg) = (frag_msg*)realloc((*f_msg), ((*line_m)+1) * sizeof(frag_msg));
						if ((*f_msg) == NULL) { fprintf(stderr, "frag dp path] Not enough memory.\n"); exit(0); }
						(*line_m)++;
					}
					frag_copy_msg((*f_msg), (*f_msg)+cur_num);
					cur_line=cur_num; frag_num = 0;
					//frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
					cur_num++;

					if (cur_num == (*line_m))
					{
						(*f_msg) = (frag_msg*)realloc((*f_msg), ((*line_m)+1) * sizeof(frag_msg));
						 if ((*f_msg) == NULL) { fprintf(stderr, "[frag dp path] Not enough memory.\n"); exit(0); }
						 (*line_m)++;
					}
					frag_copy_msg((*f_msg), (*f_msg)+cur_num);
					frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, f_node, _line, _m_len);
					cur_num++;
				}
				//else frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
			}
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
			// MIS or REF-N
			//new frag for REF-N
		}*/
		else if (f_node[cur_x][cur_y].match_flag == F_MISMATCH)
		{
			if (abs(pre_x - cur_x) > 1)
			{
				//dp-line for INV/TRS
				_m_len = 0;
				_m_len = frag_mini_dp_line(f_node, a_msg, seed_len, (line_node){pre_x, pre_y}, (line_node){cur_x, cur_y}, _line, 0, 0, n_seed);
				if (_m_len != 0)
				{
					frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
					((*f_msg)+cur_line)->frag_left_bound = a_msg[_line[_m_len-1].x].read_id;
					//++frag_num;
					//new frag for main-line
					if (cur_num == (*line_m))
					{
						(*f_msg) = (frag_msg*)realloc((*f_msg), ((*line_m)+1) * sizeof(frag_msg));
						if ((*f_msg) == NULL) { fprintf(stderr, "frag dp path] Not enough memory.\n"); exit(0); }
						(*line_m)++;
					}
					frag_copy_msg((*f_msg), (*f_msg)+cur_num);
					cur_line=cur_num; frag_num = 0;
					frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
					((*f_msg)+cur_line)->frag_right_bound = a_msg[_line[0].x].read_id;
					cur_num++;

					//new frag for INV/TRS
					if (cur_num == (*line_m))
					{
						(*f_msg) = (frag_msg*)realloc((*f_msg), ((*line_m)+1) * sizeof(frag_msg));
						 if ((*f_msg) == NULL) { fprintf(stderr, "[frag dp path] Not enough memory.\n"); exit(0); }
						 (*line_m)++;
					}
					frag_copy_msg((*f_msg), (*f_msg)+cur_num);
					frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, f_node, _line, _m_len);
					((*f_msg)+cur_num)->frag_left_bound = a_msg[pre_x].read_id;
					((*f_msg)+cur_num)->frag_right_bound = a_msg[cur_x].read_id;
					cur_num++;
				}
				//XXX
				else frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
			}
			//MIS, REF-N
			//no seeds' aln-res exist
			else
			{
				//frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
				//++frag_num;
				//frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+cur_line, frag_num, seed_len);
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
			}
		}
		else if (f_node[cur_x][cur_y].match_flag == F_MATCH)	//: MATCH
		{
			frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
		}
		else {fprintf(stderr, "[frag dp path] Error: Unknown flag, \"%d\"", f_node[cur_x][cur_y].match_flag); exit(0);}
	}
	//last start
	cur_x = line[0].x; cur_y = line[0].y;
	frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+cur_line, frag_num, seed_len);
	int left_bound = 0;
	//MIS-MATCH
	if (cur_x > 0)
	{
		_m_len = 0;
		_m_len = frag_mini_dp_line(f_node, a_msg, seed_len, (line_node){-1, 0}, (line_node){cur_x, cur_y}, _line, 0, 0, n_seed);
		if (_m_len != 0)
		{
			if (cur_num == (*line_m))
			{
				(*f_msg) = (frag_msg*)realloc(*f_msg, ((*line_m)+1) * sizeof(frag_msg));
				if ((*f_msg) == NULL) { fprintf(stderr, "[frag_dp_path] Not enough memory.(line_m: %d)\n", (*line_m)+1); exit(0); }
				(*line_m)++;
			}
			frag_copy_msg(*f_msg, (*f_msg)+cur_num);
			frag_mini_dp_path(a_msg, seed_len, (*f_msg)+cur_num, f_node, _line, _m_len);
			((*f_msg)+cur_num)->frag_left_bound = 0;
			((*f_msg)+cur_num)->frag_right_bound = a_msg[cur_x].read_id;
			cur_num++;
			left_bound = a_msg[_line[_m_len-1].x].read_id;
		}
	}
	((*f_msg)+cur_line)->frag_left_bound = left_bound;
	//free(_line);

	(*line_n) = cur_num;
	/*
	int j, k, s_i, a_i;
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
	//free variables
	//free(line); for (i = 0; i < n_seed; ++i) free(f_node[i]); free(f_node);
	return 1;
}

//@: remain frags with more than 2 seeds.
/*int backtrack(aln_msg* a_msg, path_msg **path, int n_seed, int *price_n, int seed_len, frag_msg *f_msg)  //from end to start, find every fragment's postion
{
    if (n_seed == -1)
    {
        return 0;
    }
    //Determin the start point of backtrack.
    int i;//(path[i][min_i])
    for (i = n_seed-1; i >= 0; i--)
    	if (price_n[i] > 0)	break;

    //backtrack from (i, min_i)
    int last_x = i, last_y = price_n[i]-1, frag_num = 0, tmp, flag = 0;
    int end_x, start_x;
    //first end
    frag_set_msg(a_msg, last_x, last_y, 1, f_msg, frag_num, seed_len);
    end_x = last_x;
	
	fprintf(stdout, "    end   %d ref %lld read %d\n", a_msg[last_x].at[last_y].nsrand, (long long)a_msg[last_x].at[last_y].offset, a_msg[last_x].read_id);
    while (path[last_x][last_y].flag != F_PATH_END) {
    	if (path[last_x][last_y].flag != F_MATCH) {
    		//start
    		start_x = last_x;
    		if (start_x == end_x) {
    			;fprintf(stdout, "    flase start %d ref %lld read %d\n", a_msg[last_x].at[last_y].nsrand, (long long)a_msg[last_x].at[last_y].offset, a_msg[last_x].read_id);
    		}
			else {
				frag_set_msg(a_msg, last_x, last_y, 0, f_msg, frag_num, seed_len);
				fprintf(stdout, "    start %d ref %lld read %d\n", a_msg[last_x].at[last_y].nsrand, (long long)a_msg[last_x].at[last_y].offset, a_msg[last_x].read_id);
    			++frag_num;
    		}
    		flag = 1;
    	}
    	tmp = last_x;
    	last_x = path[last_x][last_y].from.x;
    	last_y = path[tmp][last_y].from.y;

    	if (flag == 1) {
    		//next end
    		frag_set_msg(a_msg, last_x, last_y, 1, f_msg, frag_num, seed_len);
			end_x = last_x;

			fprintf(stdout, "    end   %d ref %lld read %d\n", a_msg[last_x].at[last_y].nsrand, (long long)a_msg[last_x].at[last_y].offset, a_msg[last_x].read_id);
    		flag = 0; 
    	} else {
    		frag_set_msg(a_msg, last_x, last_y, 2, f_msg, frag_num, seed_len);	//seed
    		//printf("seed: %d %d\n", a_msg[last_x].read_id, last_y);
    	}
    }
	//start
    frag_set_msg(a_msg, last_x, last_y, 0, f_msg, frag_num, seed_len);
	fprintf(stdout, "    start %d ref %lld read %d\n", a_msg[last_x].at[last_y].nsrand, (long long)a_msg[last_x].at[last_y].offset, a_msg[last_x].read_id);
	return 1;
}*/

int frag_cluster(const char *read_prefix, char *seed_result, seed_msg *s_msg, int seed_len, bntseq_t *bns, uint8_t *pac)
{
	FILE *result_p; char readline[1024];
	int n_read/*start from 1*/, n_seed, i; char srand;
	int read_id, chr, edit_dis; long long offset; char cigar[1024];
	
	aln_msg *a_msg; 
	frag_msg *f_msg; int line_n, line_m;

	gzFile readfp; kseq_t *read_seq_t; char *read_seq;

	//alloc mem and initialization
	a_msg = aln_init_msg(s_msg->seed_max);

	//f_node = (frag_DP_node*)malloc(s_msg->seed_max * sizeof(frag_DP_node));
	line_node *line = (line_node*)malloc(s_msg->seed_max * sizeof(line_node));
	line_node *_line = (line_node*)malloc(s_msg->seed_max * sizeof(line_node));
	frag_dp_node **f_node = (frag_dp_node**)malloc(s_msg->seed_max * sizeof(frag_dp_node*));
	for (i = 0; i < s_msg->seed_max; ++i)
		f_node[i] = (frag_dp_node*)malloc(PER_ALN_N * sizeof(frag_dp_node)); 
	if (line == NULL || f_node == NULL) { fprintf(stderr, "[frag_dp_path] Not enougy memory.\n"); exit(0); }
	
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
	//get seed msg of every read
	while (fgets(readline, 1024, result_p) != NULL)
	{
		//XXX for new-out.0 add 'cigar'
		sscanf(readline, "%d %d %lld %c %d %s", &read_id, &chr, &offset, &srand, &edit_dis, cigar);
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
					//if (find_path(a_msg, n_seed, line, path_end, path, price_n, seed_len, bns->n_seqs))
					//if (frag_find_path(a_msg, n_seed, seed_len, line, line_end, f_node, f_msg))
					f_msg[0].last_len = s_msg->last_len[n_read];
					f_msg[0].seed_all = s_msg->n_seed[n_read]-s_msg->n_seed[n_read-1];
					//fprintf(stdout, "%s\n", s_msg->read_name[n_read]);
					if (frag_dp_path(a_msg, n_seed, seed_len, &f_msg, &line_n, &line_m, line, f_node, _line))
					{
						//if (backtrack(a_msg, path, n_seed, price_n, seed_len, f_msg))
						/* SW-extenging */
						frag_check(s_msg->read_name[n_read], bns, pac, read_prefix, read_seq, s_msg->read_len[n_read], &f_msg, line_n, a_msg, &hash_num, &hash_node, seed_len);
					}
				}
				n_seed = 0;
				while (s_msg->n_seed[n_read] < read_id) {
					if (FLAG == 0) FLAG = 1;
					//else fprintf(stdout, "read %d n_seed 0\nfrag: 0\n\n", n_read);
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
			setAmsg(a_msg, n_seed, multi_aln, read_id-s_msg->n_seed[n_read-1], chr, offset, srand, edit_dis, cigar);
		}
	}
	//fprintf(stdout, "read %d start %d n_seed %d\n", n_read, s_msg->n_seed[n_read-1]+1, n_seed);
	//if (find_path(a_msg, n_seed, line, path_end, path, price_n, seed_len, bns->n_seqs))
	//if (frag_find_path(a_msg, n_seed, seed_len, line, line_end, f_node, f_msg))
	f_msg[0].seed_all = s_msg->n_seed[n_read] - s_msg->n_seed[n_read-1];
	f_msg[0].last_len = s_msg->last_len[n_read];
	//fprintf(stdout, "%s\n", s_msg->read_name[n_read]);
	if (frag_dp_path(a_msg, n_seed, seed_len, &f_msg, &line_n, &line_m, line, f_node, _line))
	{
		//if (backtrack(a_msg, path, n_seed, price_n, seed_len, f_msg))
		/* SW-extenging */
		frag_check(s_msg->read_name[n_read], bns, pac, read_prefix, read_seq, s_msg->read_len[n_read], &f_msg, line_n, a_msg, &hash_num, &hash_node, seed_len);
	}
	fclose(result_p);
	aln_free_msg(a_msg, s_msg->seed_max);
	for (i = 0; i < s_msg->seed_max; ++i) free(f_node[i]); free(f_node);
	free(line); free(_line);

    //hash map
	free(hash_num); 
	for (i = 0; i < hash_size; ++i) free(hash_node[i]);
	free(hash_node);

	frag_free_msg(f_msg, line_m);
	gzclose(readfp); kseq_destroy(read_seq_t);

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

int lsat_aln_core(const char *ref_prefix, const char *read_prefix, int seed_info, int no_soap2_dp, char *seed_result, char *opt_m, int opt_l)
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
	int no_soap2_dp=0, seed_info=0, opt_l=0;
	char result_f[1024]="", opt_m[100];
	char *ref, *read;
	
	opt_l = SEED_LEN;
	strcpy(opt_m, "-m 3e ");

	while ((c =getopt(argc, argv, "nsa:m:l:")) >= 0)
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

	lsat_aln_core(ref, read, seed_info, no_soap2_dp, result_f, opt_m, opt_l);
	
	free(ref); free(read);
	return 0;
}
