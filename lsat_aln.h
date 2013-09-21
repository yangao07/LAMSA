#ifndef LSAT_ALN_H
#define LSAT_ALN_H

#define READ_INIT_MAX 3000
#define SEED_INIT_MAX 1000

#define SEED_LEN 100
#define PER_ALN_N 100

#define MATCH 0
#define INSERT 1
#define DELETION 2
#define CHR_DIF 3
#define REVERSE 4
#define THRSHOLD 50

#define PRICE_DIF_CHR 3000	//相邻read比对到不同chr上的路径代价
#define PRICE_LONG_DEL 5000

#define adjest(dis) (dis>THRSHOLD?THRSHOLD:dis)


typedef struct {	//全部read包含seed数目信息
    int n_read;		//获取的read总数目			XXX contig_all -> n_read
    int *n_seed;	//存放每条read的seed数目	XXX READ_NUM -> n_seed
    int seed_max;   //contig中分割成短read的数目最大值	//XXX read_max -> seed_max
} seed_msg;

typedef struct {
	int32_t read_id;
	int32_t *chr;
	int32_t *offset;
	int8_t *nsrand;
	int8_t *edit_dis;
	int8_t n_aln;
} aln_msg;

typedef struct {
    int flag;   //MATCH INSERT DELETION
    int from;
} path_msg;

int lsat_aln(int argc, char* argv[]);

#endif
