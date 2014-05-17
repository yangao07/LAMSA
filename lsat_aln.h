#ifndef LSAT_ALN_H
#define LSAT_ALN_H

#include <stdint.h>
#define READ_INIT_MAX 3000
//#define SEED_INIT_MAX 1000

#define SEED_LEN 100
#define PER_ALN_N 100

#define HASH_MIN_LEN 1

#define EDIT_THS 3

#define FRAG_CON_STR "MPXIDCRUSE"
#define F_MATCH 0
#define F_SPLIT_MATCH 1
#define F_MISMATCH 2
#define F_INSERT 3
#define F_DELETE 4
#define F_CHR_DIF 5
#define F_REVERSE 6
#define F_UNCONNECT 7
#define F_UNMATCH 8
#define F_PATH_END 9

#define F_TRI 10
#define F_TRI_MISMATCH 11
#define F_TRI_INSERT 12
#define F_TRI_DELETE 13


#define F_SV_PEN 2
#define F_MATCH_SCORE 1

#define DP_MIN 1
#define DP_MULTI 2

#define DP_BACK_NONE 0
#define DP_BACK_PEAK 1
#define DP_BACK_LOCATED 2

#define DEL_THD 100000	//XXX
#define THRSHOLD 50

#define PRICE_DIF_CHR 100000	//相邻read比对到不同chr上的路径代价
#define PRICE_LONG_DEL 5000
#define PRICE_SKIP	20

#define FRAG_START 0
#define FRAG_SEED 1
#define FRAG_END 2

#define adjest(dis) (dis>THRSHOLD?THRSHOLD:dis)

#define SOAP2_DP_DIR "./soap2-dp"

//copy from samtools-0.1.19
#define CIGAR_STR "MIDNSHP=XB"
/*
  CIGAR operations.
 */
/*! @abstract CIGAR: M = match or mismatch*/
#define CMATCH      0
/*! @abstract CIGAR: I = insertion to the reference */
#define CINS        1
/*! @abstract CIGAR: D = deletion from the reference */
#define CDEL        2
/*! @abstract CIGAR: N = skip on the reference (e.g. spliced alignment) */
#define CREF_SKIP   3
/*! @abstract CIGAR: S = clip on the read with clipped sequence
  present in qseq */
#define CSOFT_CLIP  4
/*! @abstract CIGAR: H = clip on the read with clipped sequence trimmed off */
#define CHARD_CLIP  5
/*! @abstract CIGAR: P = padding */
#define CPAD        6
/*! @abstract CIGAR: equals = match */
#define CEQUAL      7
/*! @abstract CIGAR: X = mismatch */
#define CDIFF       8
#define CBACK       9

#define CIGAR_SHIFT	4
#define CIGAR_GEN(l,o) ((int)(l)<<CIGAR_SHIFT|(o)) 


typedef struct {	//全部read包含seed数目信息
    int n_read;		//获取的read总数目			
	char **read_name;
    int *n_seed;	//存放每条read的seed数目    index from 1
    int *read_len;	//length of read
	int *last_len;	//last_len                  index from 1
    int seed_max;   //contig中分割成短read的数目最大值	
    int read_max_len;    //max length of read
	int read_m; 
} seed_msg;

typedef struct {
	int32_t chr;
	int64_t offset;	//1-base
	int8_t nsrand;
	int8_t edit_dis;

	uint32_t *cigar;
	int cigar_len;  //default: 7 for 3-ed
	int len_dif;	//length difference between ref and read
	int cmax;
	int bmax;		//max band-width, NOT for this seed-aln, for this in-del case.
					//max of the number of inserts and the number of dels
} aln_t;
typedef struct {
	int32_t read_id;
	int8_t n_aln;
	int skip;
	aln_t *at;	
} aln_msg;

/*typedef struct {
	int same;		//
	int8_t n_aln;	//init:0
	aln_t *at;
} hash_aln_msg;*/

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

typedef struct {
	int seed_num;
	int n_line;
	int *line_i;
	int score;	//score = sum(seed_num) - sum(SV) * SV_PEN
	int from_i;
} frag_DP_node;

typedef struct {
	line_node from;
	int seed_i;
	int aln_i;
	int score;
	uint8_t match_flag;
	int dp_flag;
	int backtrack_flag;
	int node_n;
	int peak_value;
	line_node peak;
	int next_trigger_n;
	int pre_trigger_n;
    int trigger_m;
	int *trigger;	//next...pre
} frag_dp_node;

int lsat_aln(int argc, char* argv[]);

#endif
