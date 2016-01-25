#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "split_mapping.h"
#include "bntseq.h"
#include "lamsa_aln.h"
#include "frag_check.h"
#include "ksw.h"

extern char READ_NAME[1024];

	//uint64_t counter=0;
	//hash_num : number of hash-seed that have same key_int
	//hash_node : first 8*2 bit : last 8 bp, last 32-8*2 bit : pos of seed

	//hash_i: 0-2^(32-2*kmer_len)
	// 0 - 2^16
	//return the pos of kmer that matches kmer_int, and EQUAL_FLAG means EXISTS(1) or NOT(0)
	// binary search and store
int hash_search(uint32_t **hash_num, uint64_t ***hash_node, int key_int, int kmer_int, int *EQUAL_FLAG)
{
	int left=0, right=(*hash_num)[key_int]-1, middle;
	if(right == -1) {
        *EQUAL_FLAG = 0;
		return 0;
    }
	if(right == 0){
		if((((*hash_node)[key_int][0] >> 32) & 0xffff) == kmer_int) {
			*EQUAL_FLAG = 1;
			return 0;
		} else {
            *EQUAL_FLAG = 0;
            if((int)(((*hash_node)[key_int][0] >> 32) & 0xffff) > kmer_int)
                return 0;
            return 1;
        }
	}
	while(left <= right)	// bin search
	{
		middle = (left + right)>>1;
		if((((*hash_node)[key_int][middle] >> 32) & 0xffff) == kmer_int)		{
			*EQUAL_FLAG = 1;
			return middle;
		} else if((int)(((*hash_node)[key_int][middle] >> 32) & 0xffff) > kmer_int) {
			if((middle == 0) || ((int)(((*hash_node)[key_int][middle - 1] >> 32) & 0xffff) < kmer_int)) {
                *EQUAL_FLAG = 0;
				return middle;
            } else right = middle-1;
		} else left = middle + 1;
	}
    *EQUAL_FLAG = 0;
	return ((*hash_num)[key_int]);
}

//return value: if find the kmer that match key_int or not
int hash_hit(uint32_t *hash_num, uint64_t **hash_node, int *node_i, int key_int, int kmer_int)
{
	int left=0, right=hash_num[key_int]-1, middle;

	if(right == -1)		return 0;
	if(right == 0)	{
		if(((hash_node[key_int][0] >> 32) & 0xffffffff) == kmer_int) {
			*node_i = 0;
			return 1;
		} else return 0;
	}
	while(left <= right)
	{
		middle = (left + right)>>1;
		if(((hash_node[key_int][middle] >> 32) & 0xffffffff) == kmer_int) {
			*node_i = middle;
			return 1;
		} else if((int)((hash_node[key_int][middle] >> 32) & 0xffffffff) > kmer_int)
			right = middle-1;
		else left = middle + 1;
	}
	return 0;
}

int hash_calcu(int *key_int, int *kmer_int, uint8_t *seed, int hash_len, int key_len)
{
	int i;
	(*key_int) = (*kmer_int) = 0;
	for (i = 0; i < key_len; ++i) (*key_int) = (*key_int) << 2 | hash_nt4_table[(int)seed[i]];
	for (; i < hash_len; ++i) (*kmer_int) = (*kmer_int) << 2 | hash_nt4_table[(int)seed[i]];	
	return 0;
}

	//hash_node:  ----32----||--------16--------|---------16----------       
	//             kmer_int           0               pos_num
	//                           ->
	//@function:
	//hash_node:  ----32----||--------16--------|---------16----------       
	//             kmer_int    pos_start_index   pos_already_in_num
	//                           ->
	//hash_node:  ----32----||--------16--------|---------16----------       
	//             kmer_int    pos_start_index        pos_num
int init_hash_core(uint8_t *hash_seq, int hash_offset, 
				   int hash_len, int key_len, 
				   uint32_t *hash_num, uint64_t ***hash_node, int ***hash_node_num,
                   int32_t **hash_pos)
{
	int key_int, kmer_int, node_i;
	
	hash_calcu(&key_int, &kmer_int, hash_seq, hash_len, key_len);
	//node_i = hash_search(&hash_num, hash_node, key_int, kmer_int, hash_len-key_len, &EQUAL_FLAG);
	if (hash_hit(hash_num, *hash_node, &node_i, key_int, kmer_int) != 1)
    {
        int i;
        fprintf(stderr, "hash seq:\t");
        for (i = 0; i < hash_len; ++i)
            fprintf(stderr, "%c", "ACGT"[hash_seq[i]]);
        fprintf(stderr, "\n%d %d", key_int, kmer_int);
		fprintf(stderr, "\n2nd search wrong.\n"); exit(1);
    }
	(*hash_pos)[((*hash_node)[key_int][node_i] & 0xffffffff) + (*hash_node_num)[key_int][node_i]] = hash_offset;
    //printf("hash_pos: %lld %d -> %d\n", ((*hash_node)[key_int][node_i] & 0xffffffff), (*hash_node_num)[key_int][node_i], hash_offset);
    ++(*hash_node_num)[key_int][node_i];
	return 0;
}

	//@function:
    //before:
	//hash_node:  ----32----||--------16--------|---------16----------       
	//             kmer_int           0               pos_num
	//after:                          ->
	//hash_node:  ----32----||--------16--------|---------16----------       
	//             kmer_int    pos_start_index   pos_already_in_num(0)
int hash_calcu_pos_start(int hash_size, uint64_t ***hash_node, uint32_t *hash_num)
{
	int i, j;
    uint64_t tmp, pos_start=0;
	for (i = 0; i < hash_size; ++i)
	{
		for (j = 0; j < (int)hash_num[i]; ++j)
		{
			//printf("pos_start: %d %d -> %d\n", i, (*hash_node)[i][j] >> 32, pos_start);
			tmp = pos_start;
			pos_start += ((*hash_node)[i][j] & 0xffffffff);
			(*hash_node)[i][j] = ((*hash_node)[i][j] & 0xffffffff00000000) | tmp;
		}
	}
    return 0;
}

//@: calcu the number of every seed's hit positons.
int init_hash_pos_num(uint8_t *hash_seq, int hash_len, int key_len, uint32_t **hash_num, uint64_t ***hash_node)
{
	int key_int, kmer_int, node_i, EQUAL_FLAG;
	
	hash_calcu(&key_int, &kmer_int, hash_seq, hash_len, key_len);

    EQUAL_FLAG = 0;
	//binary-search and update hash table
    if ((*hash_num)[key_int] == 0) {
        node_i = 0;
    } else node_i = hash_search(hash_num, hash_node, key_int, kmer_int, &EQUAL_FLAG);

	if (EQUAL_FLAG == 1) {
    	++(*hash_node)[key_int][node_i];
    } else {
		int n = ++(*hash_num)[key_int];
        (*hash_node)[key_int] = (uint64_t*)realloc((*hash_node)[key_int], n * sizeof(uint64_t));
		if((*hash_node)[key_int] == NULL) {
			fprintf(stderr, "[init hash] ERROR: tempnode memory not enougy!\n");
			return 0;
		}
        memmove((*hash_node)[key_int]+node_i+1, (*hash_node)[key_int]+node_i, (n-1-node_i)*sizeof(uint64_t));
		(*hash_node)[key_int][node_i] = kmer_int;
		(*hash_node)[key_int][node_i] <<= 32;
		++(*hash_node)[key_int][node_i];
	}

	return 0;
}

	//hash_len = 10 (2+8)
	//key_len = 2	//fixed value
	//hash_size = 4^2 = 16	//fixed value
int init_hash(uint8_t *ref_seq, int ref_len, int hash_len, 
			  uint32_t **hash_num, uint64_t ***hash_node, int ***hash_node_num,
              int32_t **hash_pos, 
			  int key_len, int hash_size)
{
	int i;
	uint8_t *hash_seq;

	for (i = 0; i < hash_size; ++i)
		(*hash_num)[i] = 0;
	for (i = 0; i <= ref_len-hash_len; ++i)	//create hash index
	{
		//strncpy(kmer, ref_seq,	hash_len);
        hash_seq = ref_seq+i;
		init_hash_pos_num(hash_seq, hash_len, key_len, hash_num, hash_node);
	}
	hash_calcu_pos_start(hash_size, hash_node, *hash_num);
    *hash_node_num = (int**)malloc(hash_size * sizeof(int*));
    for (i = 0; i < hash_size; ++i)
        (*hash_node_num)[i] = (int*)calloc((*hash_num)[i], sizeof(int));

	for (i = 0; i <= ref_len-hash_len; ++i)
	{
		hash_seq = ref_seq+i;
		init_hash_core(hash_seq, i, hash_len, key_len, *hash_num, hash_node, hash_node_num, hash_pos);
	}
	return 0;
}

//1-mismatch betweed two adjacent nodes is allowed
//@para: bound_flag: 1 -> _head,_tail both are set as 1, 0 -> only one.
//a_i/b_i: read_i
//a_offset: ref_i - read_i
//ref_i = read_i+offset
//
//NOTE: for overlaped-DUP: ref_len = read_len+read_len-ref_len+2*hash_len
//                       
int hash_main_dis(int a_i, int a_offset, int b_i, int b_offset, lamsa_aln_para *AP, int *con_flag,
                  int ref_len, int read_len, int ref_offset)
{
    int hash_len = AP->hash_len, hash_step = AP->hash_step;
	int dis;
	if (a_i > b_i) dis = a_offset - b_offset;
	else dis = b_offset - a_offset;
	if (abs(dis) == 0) {	//match or mismatch
            if (abs(b_i - a_i) < hash_len + 2 * hash_step)	//1-mismatch seed allowed
				*con_flag = F_MATCH;
			else if (abs(b_i - a_i) < hash_len + 6 * hash_step)	//dis: 5 hash-seeds
                *con_flag = F_MISMATCH;
            else *con_flag = F_LONG_MISMATCH;
	} else { // SV or F_UNCONNECT
		if (dis > 0) *con_flag = F_DELETE;
		else if (dis >= -(abs(b_i-a_i)-hash_len)) *con_flag = F_INSERT; // unoverlaped-ins
        //else *con_flag = F_UNCONNECT;
		// overlapped
        else if (dis <= -(AP->split_len/2)) { // SV_SPLIT_LEN
            if (ref_offset > 0) { // INS region
                if (b_i>a_i) {
                    // if (read_len-ref_len+b_i+b_offset >= b_i-(a_i+hash_len-1) && read_len-(a_i+hash_len-1+a_offset)>= b_i-(a_i+hash_len-1))
                    if (read_len-ref_len+b_offset >= -(a_i+hash_len-1) && read_len-a_offset>= b_i)
                        *con_flag = F_INSERT;
                    else *con_flag = F_UNCONNECT;
                } else if (read_len-ref_len+a_offset >= -(b_i+hash_len-1) && read_len-b_offset>= a_i)
                    *con_flag = F_INSERT;
                else *con_flag = F_UNCONNECT;
            } else { // DEL/MISMATCH region
                if (b_i> a_i) {
                    // if (b_i+b_offset >= b_i-(a_i-1) && reflen-(a_i-1+a_offset) >= b_i-(a_i-1))
                    if ((b_offset >= -(a_i-1)) && (ref_len-a_offset>=b_i))
                        *con_flag = F_INSERT;
                    else *con_flag = F_UNCONNECT;
                } else {
                    if ((a_offset >= -(b_i-1)) && (ref_len-b_offset>=a_i))
                        *con_flag = F_INSERT; // overlap-allowed
                    else *con_flag = F_UNCONNECT;
                }
            }
        } else *con_flag = F_UNCONNECT;
	}
	return abs(dis);
} 

// init for dp node, head and tail node NOT include
int hash_dp_init(hash_dp_node **h_node, 
		int *hash_pos, int *start_a, int *len_a, 
		int node_i, int read_i, 
		line_node head, 
		lamsa_aln_para *AP, int dp_flag,
		int ref_len, int read_len, int ref_offset)
{
	int i;
	if (h_node[head.x][head.y].dp_flag == UNLIMITED_FLAG)
	{
		for (i = 0; i < len_a[node_i]; ++i)
		{
			h_node[node_i][i].read_i = read_i;
			h_node[node_i][i].offset = hash_pos[start_a[node_i]+i] - read_i;
			h_node[node_i][i].from = head;
			h_node[node_i][i].score = 1;//HASH_INIT_SCORE(1, 0);
			h_node[node_i][i].node_n = 1;
			h_node[node_i][i].match_flag = F_MATCH;
			h_node[node_i][i].dp_flag = dp_flag;
		}
	}
	else
	{
		int con_flag;
		for (i = 0; i < len_a[node_i]; ++i)
		{
			h_node[node_i][i].read_i = read_i;
			h_node[node_i][i].offset = hash_pos[start_a[node_i]+i] - read_i;
			hash_main_dis(h_node[head.x][head.y].read_i, h_node[head.x][head.y].offset, read_i, hash_pos[start_a[node_i]+i] - read_i, AP, &con_flag, ref_len, read_len, ref_offset);

			if (con_flag == F_UNCONNECT) {
				h_node[node_i][i].from = (line_node){-1,0};
				h_node[node_i][i].score = 0;
				h_node[node_i][i].node_n = 0;
				h_node[node_i][i].match_flag = con_flag;
				h_node[node_i][i].dp_flag = 0 - dp_flag;
			} else {
				h_node[node_i][i].from = head;
				h_node[node_i][i].score = 2 - ((con_flag <= F_MATCH_THD)?0: HASH_SV_PEN);
				h_node[node_i][i].node_n = 1;
				h_node[node_i][i].match_flag = con_flag;
				h_node[node_i][i].dp_flag = dp_flag;
			}
		}
	}
	return 0;
}

int hash_min_extend(hash_dp_node **h_node, int *len_a, int node_i, int h_len, int min_len, int dp_flag)
{
	int i, j, k;
	// (node_i, i), (j, k)
	for (i = 0; i < len_a[node_i]; ++i)
	{
		for (j = 0; j < h_len; ++j)
		{
			if (len_a[j] == min_len)
			{
				for (k = 0; k < len_a[j]; ++k)
				{
					//counter++;
					if (h_node[node_i][i].dp_flag < 0)
						continue;
					if (h_node[node_i][i].offset == h_node[j][k].offset)
					{
						h_node[node_i][i].dp_flag = dp_flag;
						goto NextCheck;
					}
				}
			}
		}
NextCheck:;
	}
	return 0;
}

//pruning XXX
int hash_dp_update(hash_dp_node **h_node, int *len_a, 
				   int node_x, int node_y, int start, lamsa_aln_para *AP, int dp_flag,
				   int ref_len, int read_len, int ref_offset)
{
	int i, j, con_flag;
	line_node max_from;
	int max_score, max_flag;

	max_from = h_node[node_x][node_y].from;
	max_score = h_node[node_x][node_y].score;
	if (h_node[node_x][node_y].dp_flag != UNLIMITED_FLAG)
	{
		for (i = node_x-1; i >= start; --i) {
			for (j = 0; j < len_a[i]; ++j) {
				if (h_node[i][j].dp_flag == dp_flag)
				{
					hash_main_dis(h_node[i][j].read_i, h_node[i][j].offset, h_node[node_x][node_y].read_i, h_node[node_x][node_y].offset, AP, &con_flag, ref_len, read_len, ref_offset);
					if (con_flag == F_UNCONNECT)
						continue;
					if (h_node[i][j].score + 1 - ((con_flag<=F_MATCH_THD)?0:HASH_SV_PEN) > max_score) {
						max_score = h_node[i][j].score + 1 - ((con_flag<=F_MATCH_THD)?0:HASH_SV_PEN);
						max_from = (line_node){i, j};
						max_flag = con_flag;
					}
				}
			}
			//当前得分是最优得分, 跳出 DP
		}
		if (max_from.x != h_node[node_x][node_y].from.x || max_from.y != h_node[node_x][node_y].from.y)
		{
			h_node[node_x][node_y].score = max_score;
			h_node[node_x][node_y].from = max_from;
			h_node[node_x][node_y].match_flag = max_flag;
			if (max_flag == F_MATCH)
				h_node[max_from.x][max_from.y].dp_flag = 0 - dp_flag;       //extend the match-node-line as long as possible, in forward direction
			h_node[node_x][node_y].node_n += h_node[max_from.x][max_from.y].node_n;
		}
	} else {
		for (i = node_x-1; i >= start; --i) {
			for (j = 0; j < len_a[i]; ++j)
			{
				if (h_node[i][j].dp_flag == dp_flag)
				{
					if (h_node[i][j].score > max_score)
					{
						max_score = h_node[i][j].score;
						max_from = (line_node){i,j};
					}
				}
			}
		}
		if (max_score > 0)
			h_node[node_x][node_y] = (hash_dp_node){max_from, -1, -1, max_score, h_node[max_from.x][max_from.y].node_n, F_MATCH, UNLIMITED_FLAG}; 
	}
	
	return 0;
}

int hash_mini_dp_init(hash_dp_node **h_node, int *len_a, 
					  int node_i, line_node head, 
					  lamsa_aln_para *AP, int mini_dp_flag,
					  int ref_len, int read_len, int ref_offset)
{
	int i;
	if (h_node[head.x][head.y].dp_flag == UNLIMITED_FLAG)
	{
		for (i = 0; i < len_a[node_i]; ++i)
		{
			h_node[node_i][i].from = head;
			h_node[node_i][i].score = 1;
			h_node[node_i][i].node_n = 1;
			h_node[node_i][i].match_flag = F_MATCH;
			h_node[node_i][i].dp_flag = mini_dp_flag;
		}
	}
	else
	{
		int con_flag;
		for (i = 0; i < len_a[node_i]; ++i)
		{
			hash_main_dis(h_node[head.x][head.y].read_i, h_node[head.x][head.y].offset, h_node[node_i][i].read_i, h_node[node_i][i].offset, AP, &con_flag, ref_len, read_len, ref_offset);
			if (con_flag == F_UNCONNECT)
			{
				h_node[node_i][i].from = (line_node){-1,0};
				h_node[node_i][i].score = 0;
				h_node[node_i][i].node_n = 0;
				h_node[node_i][i].match_flag = con_flag;
				h_node[node_i][i].dp_flag = 0 - mini_dp_flag;
			}
			else
			{
				h_node[node_i][i].from = head;
				h_node[node_i][i].score = 2 - (con_flag <= F_MATCH_THD?0 : HASH_SV_PEN);
				h_node[node_i][i].node_n = 1;
				h_node[node_i][i].match_flag = con_flag;
				h_node[node_i][i].dp_flag = mini_dp_flag;
			}
		}
	}
	return 0;
}

//line: forward
int mini_hash_main_line(hash_dp_node **h_node, 
						int *len_a, 
                        lamsa_aln_para *AP,
						line_node head, line_node tail, 
						line_node *line,
						int ref_len, int read_len, int ref_offset)
{
	int i, j, mini_dp_flag = MULTI_FLAG;
	for (i = head.x+1; i < tail.x; ++i)
		//only part of the dp-node members need to be re-inited
		hash_mini_dp_init(h_node, len_a, i, head, AP, mini_dp_flag, ref_len, read_len, ref_offset);

	h_node[tail.x][tail.y].from = head;
	h_node[tail.x][tail.y].score = 0;
	h_node[tail.x][tail.y].node_n = 0;
	h_node[tail.x][tail.y].dp_flag = mini_dp_flag;

	for (i = head.x+2; i < tail.x; ++i)
	{
		for (j = 0; j < len_a[i]; ++j)
		{
			if (h_node[i][j].dp_flag == mini_dp_flag)
				hash_dp_update(h_node, len_a, i, j, head.x+1, AP, mini_dp_flag, ref_len, read_len, ref_offset);
		}
	}
	hash_dp_update(h_node, len_a, tail.x, tail.y, head.x+1, AP, mini_dp_flag, ref_len, read_len, ref_offset);

	int node_i = h_node[tail.x][tail.y].node_n-1;
	int last_x, last_y;
	i = h_node[tail.x][tail.y].from.x;
	j = h_node[tail.x][tail.y].from.y;
	while (i != head.x)
	{
		if (node_i < 0) {fprintf(stderr, "[mini main line] bug: node_i < 0.\n");}
		line[node_i].x = i;
		line[node_i].y = j;
		--node_i;
		last_x = h_node[i][j].from.x;
		last_y = h_node[i][j].from.y;
		i = last_x;
		j = last_y;
	}
	if (node_i >= 0) {fprintf(stderr, "[mini main line] bug: node_i >= 0.\n");}
	return h_node[tail.x][tail.y].node_n;
}

//for _head and _tail both are set as 1, use the normal penalty
//for only one of _head and _tail is set as 1, use the seed_len-limit penalty
int hash_main_line(int *hash_pos, int *start_a, int *len_a, 
			       int ref_len, int read_len, int ref_offset, int hash_seed_n, 
                   lamsa_aln_para *AP,
				   hash_dp_node **h_node, line_node *line, 
				   int _head, int _tail)
{
	int i, j, node_i;
	int min_len = HASH_MIN_LEN, min_exist=0;
	line_node head={0,0}, tail={hash_seed_n+1, 0};

	//dp init
	{
		//head/tail init
		if (_head) {
			h_node[0][0] = (hash_dp_node){{-1,0}, 0-AP->hash_len, 0, 0, 0, F_MATCH, MIN_FLAG};
			if (_tail) { h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, read_len, ref_len-read_len, 0, 0, F_UNMATCH, MIN_FLAG}; } //split_offset=ref_len-read_len; }
			else h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG}; 
		} else {
			h_node[0][0] = (hash_dp_node){{-1,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG};
			if (_tail) h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, read_len, ref_len-read_len, 0, 0, F_MATCH, MIN_FLAG}; // split_offset = ref_len-read_len; } 
			else h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG};
		}

		//seed init
		for (i = 1; i <= hash_seed_n; ++i) {
			if (len_a[i] == min_len) {
				hash_dp_init(h_node, hash_pos, start_a, len_a, i, (i-1)*AP->hash_step/*read offset*/, head, AP, MIN_FLAG, ref_len, read_len, ref_offset);
				min_exist = 1;
			}
			else hash_dp_init(h_node, hash_pos, start_a, len_a, i, (i-1)*AP->hash_step/*read offset*/, head, AP, MULTI_FLAG, ref_len, read_len, ref_offset);
		}
	}
	//dp update and backtrack
	{
		//min update and backtrack
		if (min_exist) {
			//min extend
			for (i = 1; i <= hash_seed_n; ++i) {
				if (len_a[i] > min_len)
					hash_min_extend(h_node, len_a, i, hash_seed_n+2, min_len, MIN_FLAG);
			}
			//min update
			for (i = 2; i <= hash_seed_n; ++i) {
				for (j = 0; j < len_a[i]; ++j) {
					if (h_node[i][j].dp_flag == MIN_FLAG)
						hash_dp_update(h_node, len_a, i, j, 1/*update start pos*/, AP, MIN_FLAG, ref_len, read_len, ref_offset);
				}
			}
			hash_dp_update(h_node, len_a, tail.x, tail.y, 1/*update start pos*/, AP, MIN_FLAG, ref_len, read_len, ref_offset);
			//backtrack and update for remaining blank
			int mini_len;
			line_node *_line = (line_node*)malloc(hash_seed_n * sizeof(line_node));;
			line_node right, left, tmp;
			node_i = 0;

			right = tail; left = h_node[tail.x][tail.y].from;
			while (1)
			{
				if (h_node[right.x][right.y].match_flag != F_MATCH && left.x < right.x-1)// there is a seed left
				{
					mini_len = mini_hash_main_line(h_node, len_a, AP, left, right, _line, ref_len, read_len, ref_offset);
					//add mini-line to the main-line, _line is forward, but the line is reverse
					for (i = mini_len-1; i >= 0; --i)
						line[node_i++] = _line[i];
				}
				if (left.x == head.x) break;
				line[node_i++] = left;
				right = left;
				left = h_node[right.x][right.y].from;
			}

			free(_line);
			//convert line to be forward
			for (i = 0; i < node_i/2; ++i)
			{
				tmp = line[i];
				line[i] = line[node_i-1-i];
				line[node_i-1-i] = tmp;	
			}
            return node_i;
		}
		//whole-multi update
		else {
			for (i = 2; i <= hash_seed_n; ++i) {
				for (j = 0; j < len_a[i]; ++j) {
					if (h_node[i][j].dp_flag == MULTI_FLAG)
						hash_dp_update(h_node, len_a, i, j, 1/*update start pos*/, AP, MULTI_FLAG, ref_len, read_len, ref_offset);
				}
			}
			hash_dp_update(h_node, len_a, tail.x, tail.y, 1/*update start pos*/, AP, MULTI_FLAG, ref_len, read_len, ref_offset);
			node_i = h_node[tail.x][tail.y].node_n-1;

			int last_x, last_y;
			i = h_node[tail.x][tail.y].from.x;
			j = h_node[tail.x][tail.y].from.y;
			while (i != head.x)
			{
				if (node_i < 0) { fprintf(stderr, "[hash main line] bug: node_i < 0. %d %d\n", _head, _tail); exit(1); }
				line[node_i].x = i;
				line[node_i].y = j;
				--node_i;
				last_x = h_node[i][j].from.x;
				last_y = h_node[i][j].from.y;
				i = last_x;
				j = last_y;
			}
			if (node_i >= 0) {fprintf(stderr, "[hash main line] bug: node_i >= 0.\n"); exit(1);}
			return h_node[tail.x][tail.y].node_n;
		}
	}
}

//return the overlap len of two nodes 
//return value >= 0
int make_indel_cigar(int ref_left, int read_left, int ref_right, int read_right, 
                     int *clen, cigar32_t **cigar, int split_len, int *split_flag)
{
	int dlen, ilen;
	dlen = ref_left - ref_right + 1;
	ilen = read_left - read_right + 1;
	if ((dlen < 0) && (ilen < 0))
	{
		fprintf(stderr, "[make_indel_cigar] Error: dlen: %d, ilen: %d.\n", dlen, ilen); exit(1);
	}
	int len = ilen - dlen;
	if (len > 0)
	{
		(*clen) = 1;
		(*cigar)[0] = (len << 4) + CDEL;
        if (len >= split_len) *split_flag |= 2;
	}
	else if (len < 0)
	{
		(*clen) = 1;
		(*cigar)[0] = ((0-len) << 4) + CINS;
        if ((-len) >= split_len) *split_flag |= 2;
	}
	else //len==0
		(*clen) = 0;
	return (dlen > ilen ? dlen : ilen); 
}

int hash_split_map(cigar32_t **split_cigar, int *split_clen, int *split_m,
				   uint8_t *ref_seq, int ref_len, int ref_offset/*for DUP*/, uint8_t *read_seq, int read_len, 
                   lamsa_aln_para *AP,
				   uint32_t *hash_num, uint64_t **hash_node, int **hash_node_num, int32_t *hash_pos, 
				   int _head, int _tail)
{
    int hash_len = AP->hash_len, hash_step = AP->hash_step, key_len = AP->hash_key_len, split_len = AP->split_pen;
	int i, j, res=0;   //res: 0x01 -> pull trigger & cut cigar/ 0x10 -> split result
	uint8_t *read_query;
	int q_key_int, q_kmer_int, node_i;

    int hash_seed_n = (read_len-hash_len)/hash_step + 1;
	int *start_a = (int*)malloc((hash_seed_n + 2) * sizeof(int));	//2: for head and tail
	int *len_a = (int*)malloc((hash_seed_n + 2) * sizeof(int));
	if (start_a == NULL || len_a == NULL) {
        fprintf(stderr, "[hash_split_map] Not enougy memory.(ref_len %d, read_len %d)\n", ref_len, read_len); exit(1);
    }

	(*split_clen) = 0;

	//hash map, restore hash map result
	{
		len_a[0] = 1;		//for head node
		for (i = 0; i <= read_len - hash_len; i+=hash_step)
		{
			read_query = read_seq + i;
			if (hash_calcu(&q_key_int, &q_kmer_int, read_query, hash_len, key_len) == 1)	// hash 'N'
			{
				start_a[i/hash_step+1] = 0;
				len_a[i/hash_step+1] = 0;
			} else {
				if (hash_hit(hash_num, hash_node, &node_i, q_key_int, q_kmer_int) == 1)
				{
					start_a[i/hash_step + 1] = ((hash_node[q_key_int][node_i] & 0xffffffff));	//for hit seeds
                    // hash_max_hits: 50
					if ((len_a[i/hash_step + 1] = hash_node_num[q_key_int][node_i]) > 50) len_a[i/hash_step+1] = 0;
				}
				else len_a[i/hash_step + 1] = 0;	//for un-hit seeds
			}
		}
		len_a[i/hash_step + 1] = 1;	//for tail node
	}
	//hash-node DP, select and arrange seeds's order on ref
    line_node *line; hash_dp_node **h_node; int m_len;
	{
		line = (line_node*)malloc(hash_seed_n * sizeof(line_node));	//store the path of DP result
		h_node = (hash_dp_node **)malloc((hash_seed_n + 2) * sizeof(hash_dp_node*));
		for (i = 0; i < hash_seed_n+2; ++i)
			h_node[i] = (hash_dp_node*)malloc(len_a[i] * sizeof(hash_dp_node));
		if (line == NULL || h_node == NULL) {fprintf(stderr, "[hash_split_map] Not enougy memory.\n"); exit(1);} 

		m_len = hash_main_line(hash_pos, start_a, len_a, ref_len, read_len, ref_offset, hash_seed_n, AP, h_node, line, _head, _tail);
	}

	int _q_len, _t_len, _clen=0, _cm, _b_w; 
	cigar32_t *_cigar=0;
	int tail_in=hash_len/2, head_in=(hash_len+1)/2; // tail_in = (L)/2, head_in = (L+1)/2
    if (m_len > 0) {
		//fill blank with generated SV and SW, return the whole cigar
		cigar32_t *g_cigar;
		g_cigar = (cigar32_t*)malloc(sizeof(cigar32_t));
		//1. fix the region between left bound and first line
		int _refi = h_node[line[0].x][line[0].y].read_i + h_node[line[0].x][line[0].y].offset;
		int _readi = h_node[line[0].x][line[0].y].read_i;
		_q_len = _readi + tail_in;
		_t_len = _refi + tail_in;
		if (_head) {
			if (_readi != 0 && _refi != 0) { // blank exists
				if (_t_len < AP->split_len && _q_len < AP->split_len) {
					_b_w = abs(_t_len-_q_len)+3;
					ksw_global2(_q_len, read_seq, _t_len, ref_seq, 5, AP->sc_mat, AP->del_gapo, AP->del_gape, AP->ins_gapo, AP->ins_gape, _b_w, &_clen, &_cigar);
				} else {
					_b_w = (abs(_t_len-_q_len) > hash_len) ? hash_len : (abs(_t_len-_q_len)+3);
					res |= ksw_bi_extend(_q_len, read_seq, _t_len, ref_seq, 5, AP->sc_mat, _b_w, hash_len*AP->match, hash_len*AP->match, AP, &_cigar, &_clen, &_cm);
				} 					
				_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
				free(_cigar);
			} else {//no blank, add SV cigar
				//         ref_left, read_left,  ref_right, read_right
				make_indel_cigar(-1, -1, _refi, _readi, &_clen, &g_cigar, split_len, &res);
				_push_cigar(split_cigar, split_clen, split_m, g_cigar, _clen);
				_push_cigar1(split_cigar, split_clen, split_m, (tail_in<<4)|CMATCH);
			}
		}
        //2. fix the region between match-lines
		int start_i, overlap = 0;	//F_UNMATCH seeds' overlap
		start_i = 0;
		for (i = 0; i < m_len; ++i) {
			if (i == m_len-1 || h_node[line[i+1].x][line[i+1].y].match_flag >= F_MATCH_THD) {
				//start -> i
				g_cigar[0] = (h_node[line[i].x][line[i].y].read_i-h_node[line[start_i].x][line[start_i].y].read_i+hash_len-tail_in-head_in-overlap) << 4 | CMATCH;
				_push_cigar1(split_cigar, split_clen, split_m, g_cigar[0]);
				if (i == m_len-1) break;
				int l_readi = h_node[line[i].x][line[i].y].read_i + hash_len-1;
				int r_readi = h_node[line[i+1].x][line[i+1].y].read_i;
				int l_refi  = h_node[line[i].x][line[i].y].read_i + hash_len + h_node[line[i].x][line[i].y].offset - 1;
				int r_refi  = h_node[line[i+1].x][line[i+1].y].read_i + h_node[line[i+1].x][line[i+1].y].offset;
				int l_offset = h_node[line[i].x][line[i].y].offset;
				int r_offset = h_node[line[i+1].x][line[i+1].y].offset;
				if (l_readi + 1 < r_readi && l_refi + 1 < r_refi) { // blank exists
					_q_len = r_readi - (l_readi + 1) + head_in + tail_in;
					_t_len = _q_len + r_offset - l_offset;

					if (_q_len < AP->split_len && _t_len < AP->split_len) {
						_b_w = abs(_q_len-_t_len)+3;
						ksw_global2(_q_len, read_seq+l_readi+1-head_in, _t_len, ref_seq+l_refi+1-head_in, 5, AP->sc_mat, AP->del_gapo, AP->del_gape, AP->ins_gapo, AP->ins_gape, _b_w, &_clen, &_cigar);
					} else {
						_b_w = (abs(_t_len-_q_len) > hash_len) ? hash_len : (abs(_t_len-_q_len)+3);
						res |= ksw_bi_extend(_q_len, read_seq+l_readi+1-head_in, _t_len, ref_seq+l_refi+1-head_in, 5, AP->sc_mat, _b_w, hash_len*AP->match, hash_len*AP->match, AP, &_cigar, &_clen, &_cm);
					} 
					_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
					overlap = 0;
					free(_cigar);
				} else if (l_refi >= r_refi) { // overlap exists
					// similar to the codes in [split_mapping](frag_check.c)
					_q_len = r_readi - (l_readi + 1) + head_in;
					_t_len = _q_len + (ref_offset>0?hash_len:0);
					int lqe, lte, rqe, rte;
					// left-extend
					ksw_extend_core(_q_len, read_seq+l_readi+1-head_in, _t_len, ref_seq+l_refi+1-head_in, 5, AP->sc_mat, 3, hash_len*AP->match, AP, &lqe, &lte, &_cigar, &_clen, &_cm);
					_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
					free(_cigar);
                    // right-extend
                    _q_len = r_readi - (l_readi + 1) + tail_in;
                    _t_len = _q_len + (ref_offset>0?hash_len:0);
                    uint8_t *re_qseq = (uint8_t*)malloc(_q_len * sizeof(uint8_t));
                    uint8_t *re_tseq = (uint8_t*)malloc(_t_len * sizeof(uint8_t));
                    if (r_readi+tail_in-_q_len < 0 || r_refi+tail_in-_t_len < -ref_offset-(ref_offset>0?hash_len:0)) {
                        fprintf(stderr, "[hash_split_map] BUG. %s\n", READ_NAME); exit(1);
                    }
                    for (j = 0; j < _q_len; ++j) re_qseq[j] = read_seq[r_readi+tail_in-1-j];
                    for (j = 0; j < _t_len; ++j) re_tseq[j] = ref_seq[r_refi+tail_in-1-j];
                    ksw_extend_core(_q_len, re_qseq, _t_len, re_tseq, 5, AP->sc_mat, 3, hash_len*AP->match, AP, &rqe, &rte, &_cigar, &_clen, &_cm);
                    _invert_cigar(&_cigar, _clen);
                    free(re_qseq); free(re_tseq);
 				
					// merge, add overlap-flag('S/H')
					int Sn = _q_len + head_in - lqe - rqe, Hn = r_refi + head_in+tail_in - l_refi - 1 - lte - rte;
                    _push_cigar0(split_cigar, split_clen, split_m, (Sn<<4)|CSOFT_CLIP); 
                    _push_cigar0(split_cigar, split_clen, split_m, (Hn<<4)|CHARD_CLIP);
                    _push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
                    overlap = 0;
                    free(_cigar);
                } else {// no blank & overlap, push SV cigar
                    _push_cigar1(split_cigar, split_clen, split_m, (head_in<<4)|CMATCH);
                    overlap = make_indel_cigar(l_refi, l_readi, r_refi, r_readi, &_clen, &g_cigar, split_len, &res);
                    _push_cigar(split_cigar, split_clen, split_m, g_cigar, _clen);
					_push_cigar1(split_cigar, split_clen, split_m, (tail_in<<4)|CMATCH);
				}
				start_i = i+1;
			}
		}
		//3. fix the region between last line and right bound	//i = line_end[m_i]-1
		_readi = h_node[line[m_len-1].x][line[m_len-1].y].read_i+hash_len-1;
		_refi = h_node[line[m_len-1].x][line[m_len-1].y].read_i+h_node[line[m_len-1].x][line[m_len-1].y].offset + hash_len-1;
		_q_len = read_len - (_readi+1)+head_in;
		_t_len = ref_len - (_refi+1)+head_in;
		if (_tail) {
			if (_readi + 1 < read_len && _refi + 1 < ref_len) { // blank exists
				if (_q_len < AP->split_len && _t_len < AP->split_len) {
					_b_w = abs(_t_len-_q_len)+3;
					ksw_global2(_q_len, read_seq+_readi+1-head_in, _t_len, ref_seq+_refi+1-head_in, 5, AP->sc_mat, AP->del_gapo, AP->del_gape, AP->ins_gapo, AP->ins_gape, _b_w, &_clen, &_cigar);
				} else {
					_b_w = (abs(_t_len-_q_len) > hash_len) ? hash_len : (abs(_t_len-_q_len)+3);
					res |=  ksw_bi_extend(_q_len, read_seq+_readi+1-head_in, _t_len, ref_seq+_refi+1-head_in, 5, AP->sc_mat, _b_w, hash_len*AP->match, hash_len*AP->match, AP, &_cigar, &_clen, &_cm);
				}
				_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
				free(_cigar);
			} else { // no blank, add SV cigar
				_push_cigar1(split_cigar, split_clen, split_m, (head_in<<4)|CMATCH);
				make_indel_cigar(_refi, _readi, ref_len, read_len, &_clen, &g_cigar, split_len, &res);
				_push_cigar(split_cigar, split_clen, split_m, g_cigar, _clen);
			}
		} 
		free(g_cigar);
	} else { // no hash-dp line nodes exist
		if (_head && _tail) {
			_t_len = ref_len; _q_len = read_len; 
			if (_t_len < AP->split_len && _q_len < AP->split_len) {
				_b_w = abs(_t_len-_q_len)+3;
				ksw_global2(_q_len, read_seq, _t_len, ref_seq, 5, AP->sc_mat, AP->del_gapo, AP->del_gape, AP->ins_gapo, AP->ins_gape, _b_w, &_clen, &_cigar);
			} else {
				_b_w = (abs(_t_len-_q_len) > hash_len) ? hash_len : (abs(_t_len-_q_len)+3);
				res |= ksw_bi_extend(_q_len, read_seq, _t_len, ref_seq, 5, AP->sc_mat, _b_w, hash_len*AP->match, hash_len*AP->match, AP, &_cigar, &_clen, &_cm);
			}
			_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
			free(_cigar);
		}
		//_head or _tail -> left or right
	}
	for (i = 0; i < hash_seed_n+2; ++i) free(h_node[i]);
	free(h_node); free(line); free(start_a); free(len_a);
	return res;
}

//insertion length = read_len - ref_len
// return: 0x00->normal/ 0x01->trigger/ 0x10->split result
int split_indel_map(cigar32_t **res_cigar, int *res_len, int *res_m,
                    uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, 
                    int ref_offset,
                    lamsa_aln_para *AP,
                    uint32_t **hash_num, uint64_t ***hash_node)
{
    int hash_len = AP->hash_len, key_len = AP->hash_key_len, hash_size = AP->hash_size;

	int32_t *hash_pos = (int32_t*)malloc(ref_len * sizeof(int32_t)); int **hash_node_num;
	init_hash(ref_seq, ref_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, key_len, hash_size);

	int res;
	res = hash_split_map(res_cigar, res_len, res_m, ref_seq, ref_len, ref_offset, read_seq, read_len, AP, *hash_num, *hash_node, hash_node_num, hash_pos, 1, 1);

	free(hash_pos);
	int i; for (i = 0; i < hash_size; ++i) free(hash_node_num[i]);
	free(hash_node_num);
	return res;
}

int hash_right_bound_map(cigar32_t **cigar, int *cigar_len, int *cigar_m,
		uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
        lamsa_aln_para *AP,
		uint32_t **hash_num, uint64_t ***hash_node)
{
    int hash_len = AP->hash_len, hash_key = AP->hash_key_len, hash_size = AP->hash_size;
	int32_t *hash_pos = (int32_t*)malloc(ref_len * sizeof(int32_t));
	int **hash_node_num;
	if (init_hash(ref_seq, ref_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, hash_key, hash_size) != 0)
	{
		(*cigar_len) = 0;
		free(hash_pos);
		return 1;
	}
	int res = 0;
	res |= hash_split_map(cigar, cigar_len, cigar_m, ref_seq, ref_len, 0, read_seq, read_len, AP, *hash_num, *hash_node, hash_node_num, hash_pos, 1, 0);
	int i, readINcigar=0, refINcigar=0;
	for (i = 0; i < (*cigar_len); ++i)
	{
		if (((*cigar)[i] & 0xf) == CMATCH) readINcigar += ((*cigar)[i] >> 4), refINcigar += ((*cigar)[i] >> 4);
		else if (((*cigar)[i] & 0xf) == CINS || ((*cigar)[i] & 0xf) == CSOFT_CLIP) readINcigar += ((*cigar)[i] >> 4);
		else if (((*cigar)[i] & 0xf) == CDEL || ((*cigar)[i] & 0xf) == CHARD_CLIP) refINcigar += ((*cigar)[i] >> 4);
		else {
			fprintf(stderr, "[hash_right_bound_map] cigar error: "); printcigar(stderr, *cigar, *cigar_len); fprintf(stderr, "\n"); exit(1);
		}
	}
	if (readINcigar < read_len)     
	{
		int qlen = read_len-readINcigar;
		int tlen = read_len -readINcigar + ((hash_len+read_len-readINcigar)*AP->match-5-1)/2;
		tlen = tlen < (ref_len - refINcigar) ? tlen : (ref_len - refINcigar);
		int read_end, ref_end, n_cigar_, m_cigar_;
		cigar32_t *cigar_= NULL;
		//
		//printf("ref:\t");for (i = 0; i < tlen; ++i) printf("%c", "ACGT"[(ref_seq+refINcigar)[i]]);printf("\n");
		//printf("read:\t");for (i =0; i < qlen; ++i) printf("%c", "ACGT"[(read_seq+readINcigar)[i]]);printf("\n");
		res |= ksw_extend_c(qlen, read_seq+readINcigar, tlen , ref_seq+refINcigar, 5, AP->sc_mat, abs(read_len-readINcigar-ref_len + refINcigar)+3, hash_len*AP->match, AP, &read_end, &ref_end, &cigar_, &n_cigar_, &m_cigar_);
		if (cigar_ != NULL)
		{
			_push_cigar(cigar, cigar_len, cigar_m, cigar_, n_cigar_);
			//for (i = 0; i < n_cigar_; ++i) (*cigar)[(*cigar_len)++] = cigar_[i];
			//printf("tail:\t"); printcigar(cigar_, n_cigar_); printf("\n");
			free(cigar_);
		}
		if (read_end < read_len - readINcigar) 
			_push_cigar1(cigar, cigar_len, cigar_m, (((read_len - readINcigar-read_end) << 4) | CSOFT_CLIP));//'S' exist
		//(*cigar)[(*cigar_len)++] = (((read_len - readINcigar-read_end) << 4) | CSOFT_CLIP); 
	}

	free(hash_pos);
	for (i = 0; i < hash_size; ++i) free(hash_node_num[i]);
	free(hash_node_num);
	return res;
}

int hash_left_bound_map(cigar32_t **cigar, int *cigar_len, int *cigar_m,
		uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
        lamsa_aln_para *AP,
		uint32_t **hash_num, uint64_t ***hash_node)
{
    int hash_len = AP->hash_len, hash_key = AP->hash_key_len, hash_size = AP->hash_size; 
	int hash_cigar_len=0, hash_cigar_m=CIGAR_LEN_M;
	cigar32_t *hash_cigar = (cigar32_t*)malloc(hash_cigar_m * sizeof(cigar32_t));;
	int32_t *hash_pos = (int32_t*)malloc(ref_len * sizeof(int32_t));
	(*cigar_len) = 0;
	int **hash_node_num;
	if (init_hash(ref_seq, ref_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, hash_key, hash_size) != 0) {
		(*cigar_len) = 0;
		free(hash_pos); free(hash_cigar);
		return 1;
	}
	int res, i, readINcigar=0, refINcigar=0;
	res = hash_split_map(&hash_cigar, &hash_cigar_len, &hash_cigar_m, ref_seq, ref_len, 0, read_seq, read_len, AP, *hash_num, *hash_node, hash_node_num, hash_pos, 0, 1);
	for (i = 0; i < hash_cigar_len; ++i)
	{
		if ((hash_cigar[i] & 0xf) == CMATCH) readINcigar+=(hash_cigar[i] >> 4), refINcigar+=(hash_cigar[i] >> 4);
		else if ((hash_cigar[i] & 0xf) == CINS || (hash_cigar[i] & 0xf) == CSOFT_CLIP) readINcigar += (hash_cigar[i] >> 4);
		else if ((hash_cigar[i] & 0xf) == CDEL || (hash_cigar[i] & 0xf) == CHARD_CLIP) refINcigar += (hash_cigar[i] >> 4);
		else fprintf(stderr, "[hash_left_bound_map] cigar error: "), printcigar(stderr, hash_cigar, hash_cigar_len), fprintf(stderr, "\n"), exit(1);
	}
	if (readINcigar < read_len)     //'S' exists
	{
		int qlen = read_len - readINcigar;
		int tlen = read_len - readINcigar + ((hash_len+read_len-readINcigar)*AP->match-5-1)/2;
		tlen = tlen < (ref_len -refINcigar) ? tlen : (ref_len - refINcigar);
		uint8_t *qseq = (uint8_t*)malloc(qlen * sizeof(uint8_t));
		uint8_t *tseq = (uint8_t*)malloc(tlen * sizeof(uint8_t));
		int read_end, ref_end, n_cigar_, m_cigar_;
		cigar32_t *cigar_= 0;
		for (i = 0; i < qlen; ++i) qseq[i] = read_seq[qlen-i-1];
		for (i = 0; i < tlen; ++i) tseq[i] = ref_seq[ref_len-refINcigar-i-1];

		res |= ksw_extend_c(qlen, qseq, tlen, tseq, 5, AP->sc_mat, abs(qlen-tlen)+3, hash_len*AP->match, AP, &read_end, &ref_end, &cigar_, &n_cigar_, &m_cigar_);

		if (cigar_ != NULL) {
			if (read_end < read_len-readINcigar) _push_cigar1(&cigar_, &n_cigar_, &m_cigar_, (((read_len-readINcigar-read_end) << 4) | CSOFT_CLIP));//'S' exsit
			//invert cigar
			_invert_cigar(&cigar_, n_cigar_);
			_push_cigar(cigar, cigar_len, cigar_m, cigar_, n_cigar_);
			free(cigar_);
		} else
			_push_cigar1(cigar, cigar_len, cigar_m, (((read_len-readINcigar) << 4) | CSOFT_CLIP)); //'S' exsit
		free(tseq); free(qseq); 
	}
	_push_cigar(cigar, cigar_len, cigar_m, hash_cigar, hash_cigar_len);

	free(hash_pos); free(hash_cigar);
	for (i = 0; i < hash_size; ++i) free(hash_node_num[i]);
	free(hash_node_num);
	return res;
}

//for 'MIS-MATCH' case
int hash_both_bound_map(cigar32_t **cigar, int *cigar_len, int *cigar_m,
		uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
        lamsa_aln_para *AP,
		uint32_t **hash_num, uint64_t ***hash_node)
{
    int hash_len = AP->hash_len, hash_key = AP->hash_key_len, hash_size = AP->hash_size;
	int32_t *hash_pos = (int32_t*)malloc(ref_len * sizeof(int32_t));
	int **hash_node_num;
	if (init_hash(ref_seq, ref_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, hash_key, hash_size) != 0) {
		(*cigar_len = 0);
		free(hash_pos);
		return 1;
	}
	int i, res;
	res = hash_split_map(cigar, cigar_len, cigar_m, ref_seq, ref_len, 0, read_seq, read_len, AP, *hash_num, *hash_node, hash_node_num, hash_pos, 1, 1);
	free(hash_pos);
	for (i = 0; i < hash_size; ++i) free(hash_node_num[i]);
	free(hash_node_num);
	return res;
}
