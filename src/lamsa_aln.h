#ifndef LAMSA_ALN_H
#define LAMSA_ALN_H

#include <stdlib.h>
#include <stdint.h>
#include "kstring.h"

#define CHUNK_SIZE 10000000
#define CHUNK_READ_N 512

#define READ_MAX_NUM 1000

//default aln_para
#define SEED_LEN 50
#define SEED_STEP 100
#define SEED_PER_LOCI 200
#define SV_MAX_LEN 10000
#define OVLP_RAT 0.7
#define MAX_SKEL 10
#define MAX_BWT_REG 300
#define PB_MIN_BWT_REG 50
#define ON_MIN_BWT_REG 100
#define BWT_KMER 19

#define ED_RATE 0.04
#define ID_RATE 0.04
#define MIS_RATE 0.04
#define MATCH_DIS 5
#define MISMATCH_THD 10
#define MAT_RATE 0.80

#define MAT_SCORE 1
#define MIS_PEN 3
#define OPEN_PEN 5
#define EXT_PEN 2
#define EXT_BAND_W 200
#define END_BONUS 5

// PacBio(10%,15%,1:12:2) aln para
#define PB_SEED_LEN 50
#define PB_SEED_STEP 25

#define PB_ED_RATE 0.3
#define PB_ID_RATE 0.3
#define PB_MAT_RATE 0.7
#define PB_MIS_RATE 0.04
#define PB_MISMATCH_THD 10

#define PB_MAT_SCORE 1
#define PB_MIS_PEN 1
#define PB_INS_OPEN_PEN 1
#define PB_INS_EXT_PEN 1
#define PB_DEL_OPEN_PEN 1
#define PB_DEL_EXT_PEN 1
#define PB_INS_EXT_OPEN_PEN 2
#define PB_INS_EXT_EXT_PEN 1
#define PB_DEL_EXT_OPEN_PEN 2
#define PB_DEL_EXT_EXT_PEN 1
#define PB_END_BONUS 0

// Oxford Nanopore(%20,%30, 3:4:5)
#define ON_SEED_LEN 50
#define ON_SEED_STEP 25

#define ON_ED_RATE 0.4 // 0.45
#define ON_ID_RATE 0.1 
#define ON_MAT_RATE 0.6
#define ON_MIS_RATE 0.04
#define ON_MISMATCH_THD 10 //XXX

#define ON_MAT_SCORE 1
#define ON_MIS_PEN 1
#define ON_INS_OPEN_PEN 1
#define ON_INS_EXT_PEN 1
#define ON_DEL_OPEN_PEN 1
#define ON_DEL_EXT_PEN 1
#define ON_INS_EXT_OPEN_PEN 1
#define ON_INS_EXT_EXT_PEN 1
#define ON_DEL_EXT_OPEN_PEN 1
#define ON_DEL_EXT_EXT_PEN 1
#define ON_END_BONUS 0

#define RES_MAX_N 10
#define SPLIT_ALN_LEN 100
#define SPLIT_ALN_PEN 10

#define SEED_FIRST_ROUND_LOCI 2

//aln_per_para

//
#define START_NODE (line_node){-1,0}


/* seeds' connect relation */
//                    0123456789
#define FRAG_CON_STR "MPXLIDCRUS"
#define F_MATCH         0
#define F_SPLIT_MATCH   1
#define F_MISMATCH      2
#define F_MATCH_THD     2
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
//  [end+3] -> (line-score, x)

#define L_MERGB 0x0
#define L_NMERG 0x1
#define L_MERGH 0x2
#define L_INTER 0x4
#define L_DUMP  0x8

#define L_EXTRA 5                    // extra line-node variables
#define L_LB(l,n) ((l)[(n)].x)   // line left-boundary
#define L_RB(l,n) ((l)[(n)].y)   // line right-boundary
#define E_LB(l,n) ((l)[(n)+1].x) // extend left-boundary
#define E_RB(l,n) ((l)[(n)+1].y) // extend right-boundary
#define L_MF(l,n) ((l)[(n)+2].x) // line merge-flag
#define L_MH(l,n) ((l)[(n)+2].y) // line merge-head

#define L_LS(l,n) ((l)[(n)+3].x) // line score
#define L_BS(l,n) ((l)[(n)+3].y) // best score of merged-lines

#define L_NM(l,n) ((l)[(n)+4].x) // dis+NM

#define _L_LB(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]].x)   // line left-boundary
#define _L_RB(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]].y)   // line right-boundary
#define _E_LB(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]+1].x) // extend left-boundary
#define _E_RB(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]+1].y) // extend right-boundary
#define _L_MF(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]+2].x) // line merge-flag
#define _L_MH(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]+2].y) // line merge-head

#define _L_LS(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]+3].x) // line score
#define _L_BS(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]+3].y) // best score of merged-lines

#define _L_NM(l,sl,i) (((l)+(sl)[(i)<<1])[(sl)[((i)<<1)+1]+4].x) // dis+NM

#define _line_node(l, lsl, i) (l+(lsl)[(i)<<1])
#define _line_len(lsl, i) ((lsl)[((i)<<1)+1])


//  [end+2] -> (line-score, x)

// backtrack flag
#define DP_BACK_NONE 0
#define DP_BACK_PEAK 1
#define DP_BACK_LOCATED 2

#define FRAG_START 0
#define FRAG_SEED 1
#define FRAG_END 2

#define SOAP2_DP_DIR "./soap2-dp"

#define MAPQ_MAX 254
#define CIGAR_LEN_M 1000

/* CIGAR operations */
/* from samtools-0.1.19 */
#define CIGAR_STR "MIDNSHP=XB"
#define CIGAR_STR_HC "MIDNHHP=XB"
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

typedef struct {	
    int read_count;     //current read count
    int read_all;		//total number of reads
	int read_m; 
    int *seed_all;	    //count of seeds per read    index from 1
    int *read_len;	    //length of read
	int *last_len;	    //last_len                  index from 1
    int seed_max;       //max number of seeds in all reads
    int read_max_len;   //max length of read
} seed_msg;

typedef struct {
    char strand; int8_t nstrand;
    char chr[10]; int32_t nchr;
    int64_t offset;
    int NM;            // edit distance

    cigar_t *cigar;
    int len_dif, bmax; // max band-width, NOT for this seed-aln, for this in-del case.
					   // max of the number of inserts and the number of dels
} map_t;
typedef struct {
    map_t *map;
    int32_t map_n, map_m;
    int32_t seed_id;
    char *map_str;
} map_msg;

typedef struct {
	int readid;		// qname
	char nstrand;	// flag
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

/*typedef struct {
	int refid, is_rev;
	uint64_t ref_beg, ref_end;
    int beg, end;
} reg_t;*/

typedef struct {
    int8_t is_rev;
    int chr;
    int64_t ref_pos;
} reg_b;

typedef struct {
    reg_b *ref_beg, *ref_end;
    int beg_n, end_n, beg_m, end_m;
    int beg, end; // read
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
    int x, y;
    int a, b;
} qua_node;

typedef struct {
    int x,y,z;
} tri_node;

typedef struct {
    line_node n1, n2;
    uint8_t tri_flag; // 0x01->next, 0x10->per
} trig_node;

typedef struct {
    line_node *node;
    int *score, *NM;
    int min_score_thd;

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
	int seed_i, aln_i;      //ALLOC

	//tree generating and pruning
	int in_de; // in-degree
    int son_n; // number of son-nodes
    int son_max;            //ALLOC
    line_node *son;
    int max_score, max_NM; 
    line_node max_node;     

	int score;  // score: for connect
    int tol_NM; // total NM, from head to cur node.


	uint8_t match_flag;
	int dp_flag;

	int node_n;
} frag_dp_node;

typedef struct {
    // read msg
    //char read_name[1024];
    //int read_len;

    // seed msg
    int last_len;
    int seed_all;   //number of all the seeds
    int seed_out;   //number of seeds that are filtered out for next step
} lamsa_aln_per_para;    // specially for every read

typedef struct {
    int n_thread;   // for mulit-threads

    int seed_len, seed_step, seed_inv; // length, dis of adjacent seeds' start pos
    int per_aln_m;  // max number of seeds' aln-results
    int first_loci_thd; 

    int SV_len_thd;
    int ske_max;
    float ovlp_rat;

    int bwt_seed_len, bwt_max_len, bwt_min_len;

    int split_len;  // min length of gap that causes split-alignment
    int split_pen;  // split score penalty
    int res_mul_max;// max number of read's aln-results

    int hash_len, hash_key_len, hash_step, hash_size;

    uint8_t supp_soft, comm; // soft for supp; comment
    FILE *outp;

    int match_dis, mismatch_thd;
    int del_thd; int ins_thd; 

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
    // SW para
    int ins_gapo, ins_gape, del_gapo, del_gape; // gap open/extension penalty in sw-global
    int ins_ext_o, ins_ext_e, del_ext_o, del_ext_e; // gap open/extension penalty in sw-extend
    int match, mis;     // score matrix, match and mismatch
    int8_t sc_mat[25];
    int ext_band_w, end_bonus, zdrop;
    // read option
    float ed_rate, mis_rate, id_rate, mat_rate;
    int read_type;
    uint8_t aln_mode; // 0x3
                      //[high-id-rate][overlapping-seeding]
} lamsa_aln_para;        // for all reads
#define aln_mode_overlap_seed(m) ((m) & 1)
#define aln_mode_high_id_err(m) ((m) & 2)

int lamsa_aln(int argc, char* argv[]);
aln_reg *aln_init_reg(int read_len);
void aln_free_reg(aln_reg *a_reg);
int get_remain_reg(aln_reg *a_reg, aln_reg *re_reg, lamsa_aln_para *AP, int reg_min_thd, int reg_max_thd);
int line_pull_trg(line_node LR);

#endif
