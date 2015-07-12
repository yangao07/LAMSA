#ifndef LSAT_ALN_H
#define LSAT_ALN_H

#include <stdint.h>
#include "kstring.h"

#define READ_MAX_NUM 1000
//aln_para
#define PER_ALN_N 100
#define SPLIT_ALN_LEN 100
#define SPLIT_ALN_PEN 10
#define SV_MAX_LEN 10000
#define RES_MAX_N 10

//aln_per_para

//
#define START_NODE (line_node){-1,0}


/* seeds' connect relation */
//                    0123456789
#define FRAG_CON_STR "MPXLIDCRUS"
#define F_MATCH         0
#define F_SPLIT_MATCH   1
#define F_MISMATCH      2
#define F_LONG_MISMATCH 3

#define F_INSERT        4
#define F_DELETE        5

#define F_CHR_DIF       6
#define F_REVERSE       7
#define F_UNCONNECT     8
#define F_UNMATCH       9

#define F_TRI          10
#define F_TRI_MISMATCH 11
#define F_TRI_INSERT   12
#define F_TRI_DELETE   13

#define F_INIT 20

//node-line flag
//  [0 ~ end-1] -> node;
//  [end] -> (left-b, right-b)/(start,end)
//  [end+1] -> (extend-left-b, extend-right-b)
//  [end+2] -> (merge-flag, merge-head)
//    0x1 Not Merged(1), Merged(0)
//    0x2 Merged, Head(1) or Body(0)
//    0x4 Inter(1) or Not Inter(0)
//  [end+2] -> (line-score, x)

#define L_MERGB 0x0
#define L_NMERG 0x1
#define L_MERGH 0x2
#define L_INTER 0x4
#define L_DUMP  0x8

#define L_EXTRA 4                    // extra line-node variables
#define L_LB(l,e,i) (l[i][e[i]].x)   // line left-boundary
#define L_RB(l,e,i) (l[i][e[i]].y)   // line right-boundary
#define E_LB(l,e,i) (l[i][e[i]+1].x) // extend left-boundary
#define E_RB(l,e,i) (l[i][e[i]+1].y) // extend right-boundary
#define L_MF(l,e,i) (l[i][e[i]+2].x) // line merge-flag
#define L_MH(l,e,i) (l[i][e[i]+2].y) // line merge-head

#define L_LS(l,e,i) (l[i][e[i]+3].x) // line score
#define L_BS(l,e,i) (l[i][e[i]+3].y) // best score of merged-lines


//  [end+2] -> (line-score, x)

// backtrack flag
#define DP_BACK_NONE 0
#define DP_BACK_PEAK 1
#define DP_BACK_LOCATED 2

#define FRAG_START 0
#define FRAG_SEED 1
#define FRAG_END 2

#define SOAP2_DP_DIR "./soap2-dp"

/* CIGAR operations */
/* from samtools-0.1.19 */
#define CIGAR_STR "MIDNSHP=XB"
/* M: match or mismatch */
#define CMATCH      0
/* I: insertion to the reference */
#define CINS        1
/* D: deletion from the reference */
#define CDEL        2
/* N: skip on the reference (e.g. spliced alignment) */
#define CREF_SKIP   3
/* S: clip on the read with clipped sequence present in qseq */
#define CSOFT_CLIP  4
/* H: clip on the read with clipped sequence trimmed off */
#define CHARD_CLIP  5
/* P: padding */
#define CPAD        6
/* =: match */
#define CEQUAL      7
/* X: mismatch */
#define CDIFF       8
#define CBACK       9

#define CIGAR_SHIFT	4
#define CIGAR_GEN(l,o) ((int)(l)<<CIGAR_SHIFT|(o)) 

#define cigar32_t int32_t

typedef struct {
    cigar32_t *cigar;
    int cigar_m;
    int cigar_n;
} cigar_t;


typedef struct {	//全部read包含seed数目信息
    int read_all;		//获取的read总数目			
	char **read_name;
    int *seed_all;	//存放每条read的seed数目    index from 1
    int *read_len;	//length of read
    int *read_level;//level of read length
	int *last_len;	//last_len                  index from 1
    int seed_max;   //read中分割成的seed数目最大值	
    int read_max_len;    //max length of read
	int read_m; 
    int *seed_len;
    int *seed_inv;
} seed_msg;

typedef struct {
	int32_t chr;
	int64_t offset;	// 1-base
	int8_t nsrand;
	int NM;      // edit distance

	cigar32_t *cigar;
	int cigar_len;  // default: 7 for 3-ed
	int len_dif;	// length difference between ref and read
	int cmax;
	int bmax;		// max band-width, NOT for this seed-aln, for this in-del case.
					// max of the number of inserts and the number of dels
} aln_t;

typedef struct {
	int32_t read_id;
	int8_t n_aln;
	int skip;
	aln_t *at;	
} aln_msg;

typedef struct {
	int readid;		// qname
	char nsrand;	// flag
	int chr;		// tid
	int64_t offset;	// 1-base

	kstring_t *cigar_s;
    int NM;
	//char *cigar;	// cigar
	//int cigar_cl;	// cigar len
	//int cigar_cm;	// cigar m
} sam_t;

typedef struct {
	sam_t *sam;
	int sam_n;
	int sam_m;
} sam_msg;

typedef struct {
	int32_t x;
	int32_t y;
} from_t;

typedef struct {
    int flag;   //MATCH INSERT DELETION
    from_t from;
} path_msg;

typedef struct {
	int x;		//seed#
	int y;		//n_aln#
} line_node;
#define nodeEq(a,b) (a.x==b.x && a.y==b.y)

typedef struct {
	int refid, is_rev;
	uint64_t ref_beg, ref_end;
    int beg, end;
} reg_t;

/*typedef struct {
    // ref_beg
	int b_refid, b_is_rev;
	uint64_t ref_beg;
	// ref_end   
	int e_refid, e_is_rev;
	uint64_t ref_end;
	// read
    int beg, end;
} reg_t;*/

typedef struct {
    reg_t *reg;
    int reg_n, reg_m;
    int read_len;
} aln_reg;

typedef struct {
    int x,y,z;
} tri_node;

typedef struct {
    line_node n1, n2;
    uint8_t tri_flag; // 0x01->next, 0x10->per
} trig_node;
typedef struct {
    line_node *node;
    int *score;

    int max_n;
    int node_n;
} node_score;

#define TREE_START 0
#define TREE_BRANCH 1

typedef struct {
	int seed_num;
	int n_line;
	int *line_i;
	int score;	//score = sum(seed_num) - sum(SV) * SV_PEN
	int from_i;
} frag_DP_node;

//ALLOC //DP
typedef struct {
    int son_flag;
	line_node from; 
	int seed_i;     //ALLOC
	int aln_i;      //ALLOC

	//tree generating and pruning
	int in_de; // 入度      
    int son_n; // 子节点数目
    int son_max;            //ALLOC
    line_node *son;
    int max_score;          
    line_node max_node;     

	int score; int dis_pen; // score: for connect, dis_pen: for SV/indel penalty

	uint8_t match_flag;
	int dp_flag;

	int node_n;
    uint8_t trg_n; //00 01 10 11
    line_node next_trg, pre_trg;
} frag_dp_node;

typedef struct {
    // read msg
    char read_name[1024];
    int read_len;
    int read_level; //level of read length

    // seed msg
    int seed_len;   //length of seed
    int seed_inv;   //interval between adjacent seeds
    int seed_step;  //sum of len and inv
    int last_len;

    int seed_all;   //number of all the seeds
    int per_aln_n;  //number of seeds' aln-results
    int seed_out;   //number of seeds that are filtered out for next step

    // aln para
    int min_thd;    //minimum number of seed's aln-result for first round's DP (min-len)
    int *frag_score_table;
        /*frag_score_table[10] = {
              1,    F_MATCH
              1,    F_SPLIT_MATCH
              1,    F_MISMATCH
             -1,    F_LONG_MISMATCH
             -3,    F_INSERT
             -3,    F_DELETE
             -3,    F_CHR_DIF
             -3,    F_REVERSE
             -6,    F_UNCONNECT
             -6,    F_UNMATCH
        };*/
    int match_dis; int del_thd; int ins_thd;
} lsat_aln_per_para;    // specially for every read
//APP

typedef struct {
    int aln_type;

    int per_aln_m;  // max number of seeds' aln-results
    int SV_len_thd;
    int split_len;  // min length of gap that causes split-alignment
    int split_pen;  // split score penalty
    int res_mul_max;// max number of read's aln-results

    int hash_len, hash_key_len, hash_step, hash_size;
    int bwt_seed_len;

    // SW para
    int gapo, gape; // gap open penalty, gap extension penalty
    int match, mis;     // score matrix, match and mismatch
    int end_bonus, zdrop;
} lsat_aln_para;        // for all reads
//AP

int lsat_aln(int argc, char* argv[]);
aln_reg *aln_init_reg(int read_len);
void aln_free_reg(aln_reg *a_reg);
int get_remain_reg(aln_reg *a_reg, aln_reg *re_reg, lsat_aln_para AP, int reg_thd);
int line_pull_trg(line_node LR);

#endif
