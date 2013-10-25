#ifndef LSAT_ALN_H
#define LSAT_ALN_H

#include <stdint.h>
#define READ_INIT_MAX 3000
//#define SEED_INIT_MAX 1000

#define SEED_LEN 100
#define PER_ALN_N 100

#define EDIT_THS 3

#define MATCH 0
#define INSERT 1
#define DELETION 2
#define CHR_DIF 3
#define REVERSE 4
#define THRSHOLD 50

#define PRICE_DIF_CHR 3000	//相邻read比对到不同chr上的路径代价
#define PRICE_LONG_DEL 5000
#define PRICE_SKIP	20

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
#define CIGAR_GEN(l,o) ((l)<<CIGAR_SHIFT|(o)) 


typedef struct {	//全部read包含seed数目信息
    int n_read;		//获取的read总数目			
    int *n_seed;	//存放每条read的seed数目    index from 1
	int *last_len;	//last_len                  index from 1
    int seed_max;   //contig中分割成短read的数目最大值	
    int read_max_len;    //max length of read
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
	aln_t *at;	
} aln_msg;

typedef struct {
    int flag;   //MATCH INSERT DELETION
    int from;
	int pre_pre;	//skip one node
} path_msg;

int lsat_aln(int argc, char* argv[]);

#endif
