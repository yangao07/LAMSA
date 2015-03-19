#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "split_mapping.h"
#include "bntseq.h"
#include "lsat_aln.h"
#include "frag_check.h"
#include "ksw.h"

	//uint64_t counter=0;
	//hash_num : number of hash-seed that have same key_int
	//hash_node : first 8*2 bit : last 8 bp, last 32-8*2 bit : pos of seed

	//hash_i: 0-2^(32-2*kmer_len)
	// 0 - 2^16
	//return the pos of kmer that matches kmer_int, and EQUAL_FLAG means EXISTS(1) or NOT(0)
	// XXX binary search and store OR store and sort ?
int hash_search(uint32_t **hash_num, uint64_t ***hash_node, int key_int, int kmer_int, int kmer_len, int *EQUAL_FLAG)
{
	int left=0, right=(*hash_num)[key_int]-1, middle;
	if(right == -1)		//已存在的k-mer数量为0
    {
        *EQUAL_FLAG = 0;
		return 0;
    }
	if(right == 0)		//只存在一个k-mer
	{
		if((((*hash_node)[key_int][0] >> 32) & 0xffff) == kmer_int)
		{
			*EQUAL_FLAG = 1;
			return 0;
		}
		else
        {
            *EQUAL_FLAG = 0;
            if((((*hash_node)[key_int][0] >> 32) & 0xffff) > kmer_int)
                return 0;
            return 1;
        }
	}
	while(left <= right)	//二分查找
	{
		middle = (left + right)>>1;
		if((((*hash_node)[key_int][middle] >> 32) & 0xffff) == kmer_int)		//找到与之相等的k-mer
		{
			*EQUAL_FLAG = 1;
			return middle;
		}
		else if((((*hash_node)[key_int][middle] >> 32) & 0xffff) > kmer_int)		//大于
		{
			if((middle == 0) || ((((*hash_node)[key_int][middle - 1] >> 32) & 0xffff) < kmer_int))	//第一个大于该k-mer的k-mer
            {
                *EQUAL_FLAG = 0;
				return middle;
            }
			else
				right = middle-1;//继续二分查找
		}
		else		//小于
			left = middle + 1;
	}
    *EQUAL_FLAG = 0;
	return ((*hash_num)[key_int]);	//不存在比该k-mer大的k-mer，返回最后一个k-mer的下一位，实际上，返回值就是已存在的k-mer数量
}

	//return value: if find the kmer that match key_int or not
	//hit index XXX 
	// multi-hit XXX
int hash_hit(uint32_t *hash_num, uint64_t **hash_node, int *node_i, int key_int, int kmer_int, int kmer_len)
{
	int left=0, right=hash_num[key_int]-1, middle;

	if(right == -1)		//已存在的k-mer数量为0
		return 0;
	if(right == 0)		//只存在一个k-mer
	{
		if(((hash_node[key_int][0] >> 32) & 0xffffffff) == kmer_int)
		{
			//find the nearest hit-index
			//*hit_index = LAST_SHIFT(hash_node[key_int][0], kmer_len*2);
			//start = (hash_node[key_int][0] & 0xffff0000)>> 16;
			//len = (hash_node[key_int][0] & 0xffff);
			//*hit_index = hash_pos[start];
			//printf("hit_index: %d, len: %d, start: %d\n", *hit_index, len, start);
			/*for (i = 1; i < len; ++i)
			{
				if (abs(*hit_index-offset) > abs(hash_pos[start+i]-offset))
					*hit_index = hash_pos[start+i];
			}*/
			*node_i = 0;
			return 1;
		}
		else return 0;
	}
	while(left <= right)	//二分查找
	{
		middle = (left + right)>>1;
		if(((hash_node[key_int][middle] >> 32) & 0xffffffff) == kmer_int)		//找到与之相等的k-mer
		{
			//find the nearest hit-index
			//*hit_index = LAST_SHIFT(hash_node[key_int][middle], kmer_len*2);
			/*start = (hash_node[key_int][middle] & 0xffff0000)>> 16;
			len = (hash_node[key_int][middle] & 0xffff);
			*hit_index = hash_pos[start];
			//printf("hit_index: %d, len: %d, start: %d\n", *hit_index, len, start);
			for (i = 1; i < len; ++i)
			{
				if (abs(*hit_index-offset) > abs(hash_pos[start+i]-offset))
					*hit_index = hash_pos[start+i];
			}*/
			*node_i = middle;
			return 1;
		}
		else if(((hash_node[key_int][middle] >> 32) & 0xffffffff) > kmer_int)		//大于
			right = middle-1;//继续二分查找	
		else		//小于
			left = middle + 1;
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
				   int hash_len, int key_len, int hash_size, 
				   uint32_t *hash_num, uint64_t ***hash_node, int ***hash_node_num,
                   int32_t **hash_pos)
{
	int key_int, kmer_int, node_i;
	
	hash_calcu(&key_int, &kmer_int, hash_seq, hash_len, key_len);
	//node_i = hash_search(&hash_num, hash_node, key_int, kmer_int, hash_len-key_len, &EQUAL_FLAG);
	if (hash_hit(hash_num, *hash_node, &node_i, key_int, kmer_int, hash_len-key_len) != 1)
    {
        int i;
        fprintf(stderr, "hash seq:\t");
        for (i = 0; i < hash_len; ++i)
            fprintf(stderr, "%c", "ACGT"[hash_seq[i]]);
        fprintf(stderr, "\n%d %d", key_int, kmer_int);
		fprintf(stderr, "\n2nd search wrong.\n");
        exit(0);
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
		for (j = 0; j < hash_num[i]; ++j)
		{
			//printf("pos_start: %d %d -> %d\n", i, (*hash_node)[i][j] >> 32, pos_start);
			tmp = pos_start;
			pos_start += ((*hash_node)[i][j] & 0xffffffff);
			(*hash_node)[i][j] = ((*hash_node)[i][j] & 0xffffffff00000000) | tmp;
		}
	}
    /*printf("pos_start: %d -> %d\n", i, pos_start);
    for (i = 0; i < hash_size; ++i)
	{
		for (j = 0; j < hash_num[i]; ++j)
		{
			printf("pos_start: %d %d -> %d\n", i, (*hash_node)[i][j] >> 32, ((*hash_node)[i][j] & 0xffffffff));
			//(*hash_node)[i][j] = ((*hash_node)[i][j]&0xffffffff00000000) | (tmp<<16);
		}
	}*/
    return 0;
}

//@: calcu the number of every seed's hit positons.
int init_hash_pos_num(uint8_t *hash_seq, int hash_len, int key_len, uint32_t **hash_num, uint64_t ***hash_node)
{
	int key_int, kmer_int, node_i, EQUAL_FLAG;
	
	hash_calcu(&key_int, &kmer_int, hash_seq, hash_len, key_len);

    EQUAL_FLAG = 0;
	//binary-search and update hash table
    if ((*hash_num)[key_int] == 0)
    {
        node_i = 0;
    }
    else node_i = hash_search(hash_num, hash_node, key_int, kmer_int, hash_len-key_len, &EQUAL_FLAG);

	if (EQUAL_FLAG == 1)
    {
    	//printf("equal\n");
		//(*ha_msg)[seed_i].same = LAST_SHIFT((*hash_node)[key_int][hash_i], 2*(hash_len-key_len));
    	++(*hash_node)[key_int][node_i];
    }
    else
	{
		//printf("++\n");
		++(*hash_num)[key_int];
		//XXX
		uint64_t *tempnode = (uint64_t *)calloc((*hash_num)[key_int], sizeof(uint64_t));
		if(tempnode == NULL)
		{
			fprintf(stderr, "[init hash] ERROR: tempnode memory not enougy!\n");
			return 0;
		}
		memcpy((void*)tempnode, (const void*)(*hash_node)[key_int], (unsigned long)(node_i * sizeof(uint64_t)));		//将比新加入的k-mer小的k-mer复制至新申请的节点中

		if(tempnode == NULL)
		{
			fprintf(stderr, "[init hash] ERROR: tmpnode memory not enougy!\n");
			return 0;
		}

		tempnode[node_i] = kmer_int;		//将新加入的k-mer存入新节点
		tempnode[node_i] <<= 32;
		++tempnode[node_i];
		//tempnode[hash_i] += seed_i;
		memcpy(tempnode+node_i+1, (*hash_node)[key_int]+node_i, ((*hash_num)[key_int]-1-node_i) * sizeof(uint64_t));	//将比新加入的k-mer大的k-mer复制至新申请的节点中

		if(tempnode == NULL)
		{
			fprintf(stderr, "[hash] ERROR: tmpnode memory not enougy!\n");
			return 0;
		}

		free((*hash_node)[key_int]);		//释放旧的节点
		(*hash_node)[key_int] = tempnode;	//将新申请的节点放入hash表中
	}

	return 0;
}

	//hash_len = 10 (2+8)
	//key_len = 2	//fixed value
	//hash_size = 4^2 = 16	//fixed value
int init_hash(uint8_t *read_seq, int read_len, int hash_len, 
			  uint32_t **hash_num, uint64_t ***hash_node, int ***hash_node_num,
              int32_t **hash_pos, 
			  int key_len, int hash_size)
{
	int i;
	uint8_t *hash_seq;

	for (i = 0; i < hash_size; ++i)
		(*hash_num)[i] = 0;
	for (i = 0; i <= read_len-hash_len; ++i)	//create hash index
	{
		//strncpy(kmer, read_seq,	hash_len);
        hash_seq = read_seq+i;
		init_hash_pos_num(hash_seq, hash_len, key_len, hash_num, hash_node);
	}
	hash_calcu_pos_start(hash_size, hash_node, *hash_num);
    *hash_node_num = (int**)malloc(hash_size * sizeof(int*));
    for (i = 0; i < hash_size; ++i)
        (*hash_node_num)[i] = (int*)calloc((*hash_num)[i], sizeof(int));

	for (i = 0; i <= read_len-hash_len; ++i)
	{
		hash_seq = read_seq+i;
		init_hash_core(hash_seq, i, hash_len, key_len, hash_size, *hash_num, hash_node, hash_node_num, hash_pos);
	}
	return 0;
}

//1-mismatch betweed two adjacent nodes is allowed
//@para: bound_flag: 1 -> _head,_tail both are set as 1, 0 -> only one.
//对于每一个节点，依据于已经存在的 split_offset，设定其对应的 split_offset	XXX
int hash_main_dis(int a_i, int a_offset, int b_i, int b_offset, int hash_len, int hash_step, 
				  int *con_flag, int split_offset, 	//XXX
				  int bound_flag, int bound_len)
{
	int dis;
	if (a_i > b_i) dis = a_offset - b_offset;
	else dis = b_offset - a_offset;
	if (dis == 0)
	{	//match or mismatch
		//if (bound_flag)
		{
            if (abs(b_i - a_i) < hash_len + 2 * hash_step)	//1-mismatch seed allowed
			//if (abs(b_i - a_i) <= hash_step)	//XXX hash_step
				*con_flag = F_MATCH;
			else *con_flag = F_MISMATCH;
		}
	}
	else if (dis == split_offset)
	{
		if (dis <= (0-abs(a_i-b_i)))	//overlap-allowed
			*con_flag = F_UNCONNECT;
		else if (abs(b_i - a_i) < hash_len + 2 * hash_step)	//1-mismatch seed allowed
			*con_flag = F_SPLIT_MATCH;
		else *con_flag = F_MISMATCH;
	}
	else	//SV or F_UNCONNECT
	{
		//if (bound_flag)
		{
			if (dis > 0) *con_flag = F_INSERT;
			else if (dis < 0 && dis > (0-abs(a_i-b_i)))	//overlap-allowed
				*con_flag = F_DELETE;
			else *con_flag = F_UNCONNECT;
		}
		/*else	//XXX
		{
			if (dis > 0 && dis <= bound_len/2) *con_flag = F_INSERT;
			else if (dis < 0 && dis > (0-abs(a_i-b_i)) && dis >= 0-bound_len/2) *con_flag = F_DELETE;
			else *con_flag = F_UNCONNECT;
		}*/
		/*{
			//XXX 3 is determined by hash_len and hash_step and SW-Matrix
			if (dis > 0 && dis <= 3) *con_flag = F_INSERT;
			else if (dis < 0 && dis >= (0-3) && dis > (0-abs(a_i=b_i))) *con_flag = F_DELETE;
			else *con_flag = F_UNCONNECT;
		}*/
	}
	return abs(dis);
} 

//return the overlap len of two nodes 
//return value >= 0
int make_indel_cigar(int ref_left, int read_left, int ref_right, int read_right, int *clen, cigar32_t **cigar, int split_len, int *split_flag)
{
	int dlen, ilen;
	dlen = ref_left - ref_right + 1;
	ilen = read_left - read_right + 1;
	if ((dlen < 0) && (ilen < 0))
	{
		fprintf(stderr, "[make_indel_cigar] Error: dlen: %d, ilen: %d.\n", dlen, ilen);
		exit(-1);
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

// init for dp node, head and tail node NOT include
int hash_dp_init(hash_dp_node **h_node, 
				 int *hash_pos, int *start_a, int *len_a, 
				 int node_i, int ref_i, 
				 line_node head, 
				 int hash_len, int hash_step, int dp_flag, 
				 int split_offset, int bound_flag, int bound_len)
{
	int i;
	if (h_node[head.x][head.y].dp_flag == UNLIMITED_FLAG)
	{
		for (i = 0; i < len_a[node_i]; ++i)
		{
			h_node[node_i][i].ref_i = ref_i;
			h_node[node_i][i].offset = hash_pos[start_a[node_i]+i] - ref_i;
			h_node[node_i][i].from = head;
			//h_node[node_i][i].score = 2 * hash_len;//XXX
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
			h_node[node_i][i].ref_i = ref_i;
			h_node[node_i][i].offset = hash_pos[start_a[node_i]+i] - ref_i;
			hash_main_dis(h_node[head.x][head.y].ref_i, h_node[head.x][head.y].offset, ref_i, hash_pos[start_a[node_i]+i] - ref_i, hash_len, hash_step, &con_flag, split_offset-h_node[head.x][head.y].offset, bound_flag, bound_len);

			if (con_flag == F_UNCONNECT)
			{
				h_node[node_i][i].from = (line_node){-1,0};
				h_node[node_i][i].score = 0;
				h_node[node_i][i].node_n = 0;
				h_node[node_i][i].match_flag = con_flag;
				h_node[node_i][i].dp_flag = 0 - dp_flag;
			}
			else
			{
				h_node[node_i][i].from = head;
				//h_node[node_i][i].score = 2 * hash_len - (con_flag == F_MATCH?0: HASH_SV_PEN);//XXX
				h_node[node_i][i].score = 2 - ((con_flag<=F_SPLIT_MATCH)?0:HASH_SV_PEN);//HASH_INIT_SCORE(2, ((con_flag==F_MATCH)?0:HASH_SV_PEN));
				h_node[node_i][i].node_n = 1;
				h_node[node_i][i].match_flag = con_flag;
				h_node[node_i][i].dp_flag = dp_flag;
			}
		}
	}
	return 0;
}

//XXX use con_flag
//                                                                 h_len: whole number of hash-dp-node

#ifndef __NEW__
int hash_min_extend(hash_dp_node **h_node, int *len_a, int node_i, int h_len, int min_len, int dp_flag, int hash_len, int hash_step, int split_offset)
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
					if (h_node[node_i][i].offset == h_node[j][k].offset)    //XXX: 
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
#endif

//(node_x,node_y): the min-len node to be extended
#ifdef __NEW__
int hash_min_extend(hash_dp_node **h_node, int *len_a, 
					int node_x, int node_y, 
					int h_len, int min_len, 
					int dp_flag, 
					int hash_len, int hash_step,
					int split_offset)
{
    int i, j, con_flag;
    int last_x = node_x, last_y = node_y;

    //from right to left
    //for (i = node_x-1; i > 0; --i)
    i = node_x - 1;
    while (i > 0)
    {
        if (len_a[i] > min_len)
        {
            for (j = 0; j < len_a[i]; ++j)
            {
                //counter++;
                if (h_node[i][j].dp_flag < 0)
                    continue;
                hash_main_dis(h_node[i][j].ref_i, h_node[i][j].offset, h_node[last_x][last_y].ref_i, h_node[last_x][last_y].offset, hash_len, hash_step, &con_flag, split_offset-h_node[i][j].offset, 0, 0);
                if (con_flag == F_MATCH)
                //if (h_node[last_x][last_y].ref_i - h_node[i][j].ref_i <= hash_len+1 && h_node[i][j].offset == h_node[last_x][last_y].offset)
                {
                    h_node[i][j].dp_flag = dp_flag;
                    last_x = i;
                    last_y = j;
                    break;
                }
                else
                {
                    if (h_node[last_x][last_y].ref_i - h_node[i][j].ref_i > hash_len + hash_step)	//XXX
                        return 0;
                }
            }
        }
        --i;
    }
    last_x = node_x; last_y = node_y;
    i = node_x + 1;
    //for (i = node_x+1; i < h_len-1; ++i)
    while (i < h_len-1)
    {
        if (len_a[i] > min_len);
        {
            for (j = 0; j < len_a[i]; ++j)
            {
                if (h_node[i][j].dp_flag < 0)
                    continue;
                hash_main_dis(h_node[last_x][last_y].ref_i, h_node[last_x][last_y].offset, h_node[i][j].ref_i, h_node[i][j].offset, hash_len, hash_step, &con_flag, split_offset-h_node[last_x][last_y].offset, 0, 0);
                if (con_flag == F_MATCH)
                //if (h_node[i][j].ref_i - h_node[last_x][last_y].ref_i <= hash_len + 1 && h_node[i][j].offset == h_node[last_x][last_y].offset)
                {
                    h_node[i][j].dp_flag = dp_flag;
                    last_x = i;
                    last_y = j;
                    break;
                }
                else
                {
                    if (h_node[i][j].ref_i - h_node[last_x][last_y].ref_i > hash_len+hash_step)
                        return 0;
                }
            }
        }
        ++i;
    } 
    return 0;
}
#endif

//pruning XXX
int hash_dp_update(hash_dp_node **h_node, int *len_a, 
				   int node_x, int node_y, 
				   int start, int hash_len, int hash_step,
				   int dp_flag,
				   int split_offset, int bound_flag, int bound_len)
{
	int i, j, con_flag;
	line_node max_from;
	int max_score;
	//int score_plus, ref_plus, read_plus, pen;
	int max_flag;

	max_from = h_node[node_x][node_y].from;
	max_score = h_node[node_x][node_y].score;
	if (h_node[node_x][node_y].dp_flag != UNLIMITED_FLAG)
	{
		for (i = node_x-1; i >= start; --i)
		{
			for (j = 0; j < len_a[i]; ++j)
			{
				if (h_node[i][j].dp_flag == dp_flag)
				{
					hash_main_dis(h_node[i][j].ref_i, h_node[i][j].offset, h_node[node_x][node_y].ref_i, h_node[node_x][node_y].offset, hash_len, hash_step, &con_flag, split_offset-h_node[i][j].offset, bound_flag, bound_len);
					if (con_flag == F_UNCONNECT)
						continue;
					//overlap of pre-node and cur-node may exist
					//ref_plus = ((h_node[node_x][node_y].ref_i - h_node[i][j].ref_i) > hash_len) ? hash_len : (h_node[node_x][node_y].ref_i - h_node[i][j].ref_i);
					//read_plus = ((h_node[node_x][node_y].ref_i + h_node[node_x][node_y].offset - h_node[i][j].ref_i - h_node[i][j].offset) > hash_len) ? hash_len : (h_node[node_x][node_y].ref_i + h_node[node_x][node_y].offset - h_node[i][j].ref_i - h_node[i][j].offset);
					//pen = (con_flag == F_MATCH)?0:HASH_SV_PEN;
					//score_plus = ref_plus + read_plus - pen;
					if (h_node[i][j].score + 1 - ((con_flag<=F_SPLIT_MATCH)?0:HASH_SV_PEN) > max_score)
					//if ((HASH_DP_SCORE(h_node[i][j].score, 1, (con_flag==F_MATCH)?0:HASH_SV_PEN)) > max_score)
					{
						max_score = h_node[i][j].score + 1 - ((con_flag<=F_SPLIT_MATCH)?0:HASH_SV_PEN);
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
			if (max_flag <= F_SPLIT_MATCH)
				h_node[max_from.x][max_from.y].dp_flag = 0 - dp_flag;       //extend the mathc-node-line as long as possible, in forward direction
			h_node[node_x][node_y].node_n += h_node[max_from.x][max_from.y].node_n;
		}
	}
	else
	{
		for (i = node_x-1; i >= start; --i)
		{
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
		{
			h_node[node_x][node_y] = (hash_dp_node){max_from, -1, -1, max_score, h_node[max_from.x][max_from.y].node_n, F_MATCH, UNLIMITED_FLAG}; 
		}	
	}
	
	return 0;
}

int hash_mini_dp_init(hash_dp_node **h_node, int *len_a, 
					  int node_i, line_node head, 
					  int hash_len, int hash_step, int mini_dp_flag, 
					  int split_offset, int bound_flag, int bound_len)
{
	int i;
	if (h_node[head.x][head.y].dp_flag == UNLIMITED_FLAG)
	{
		for (i = 0; i < len_a[node_i]; ++i)
		{
			h_node[node_i][i].from = head;
			h_node[node_i][i].score = 1;	//XXX
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
			hash_main_dis(h_node[head.x][head.y].ref_i, h_node[head.x][head.y].offset, h_node[node_i][i].ref_i, h_node[node_i][i].offset, hash_len, hash_step, &con_flag, split_offset-h_node[head.x][head.y].offset, bound_flag, bound_len);
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
				h_node[node_i][i].score = 2 - (con_flag <= F_SPLIT_MATCH?0 : HASH_SV_PEN);	//XXX
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
						int *hash_pos, int *start_a, int *len_a, 
					    int hash_len, int hash_step, 
						line_node head, line_node tail, 
						line_node *line, 
						int split_offset, int bound_flag, int bound_len)
{
	int i, j, mini_dp_flag = MULTI_FLAG;
	for (i = head.x+1; i < tail.x; ++i)
		//only part of the dp-node members need to be re-inited
		hash_mini_dp_init(h_node, len_a, i, head, hash_len, hash_step, mini_dp_flag, split_offset, 1/*bound_flag*/, bound_len);

	h_node[tail.x][tail.y].from = head;
	h_node[tail.x][tail.y].score = 0;
	h_node[tail.x][tail.y].node_n = 0;
	h_node[tail.x][tail.y].dp_flag = mini_dp_flag;

	for (i = head.x+2; i < tail.x; ++i)
	{
		for (j = 0; j < len_a[i]; ++j)
		{
			if (h_node[i][j].dp_flag == mini_dp_flag)
				hash_dp_update(h_node, len_a, i, j, head.x+1, hash_len, hash_step, mini_dp_flag, split_offset, bound_flag, bound_len);
		}
	}
	hash_dp_update(h_node, len_a, tail.x, tail.y, head.x+1, hash_len, hash_step, mini_dp_flag, split_offset, bound_flag, bound_len);

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
			       int ref_len, int read_len, int hash_seed_n, 
				   int hash_len, int hash_step, 
				   hash_dp_node **h_node, line_node *line, 
				   int _head, int _tail)
{
	int i, j, node_i;
	int min_len = HASH_MIN_LEN, min_exist=0;
	line_node head={0,0}, tail={hash_seed_n+1, 0};
	int bound_flag = 0, bound_len = read_len;	//XXX
	int split_offset=read_len-ref_len;

	//dp init
	{
		//head/tail init
		if (_head) {
			h_node[0][0] = (hash_dp_node){{-1,0}, 0-hash_len, 0, 0, 0, F_MATCH, MIN_FLAG};
			if (_tail) { h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, ref_len, read_len-ref_len, 0, 0, F_UNMATCH, MIN_FLAG}; bound_flag=1;} //split_offset=read_len-ref_len; }
			else h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG}; 
		} else {
			h_node[0][0] = (hash_dp_node){{-1,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG};
			if (_tail) h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, ref_len, read_len-ref_len, 0, 0, F_MATCH, MIN_FLAG}; // split_offset = read_len-ref_len; } 
			else h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG};
		}

		//seed init
		for (i = 1; i <= hash_seed_n; ++i) {
			if (len_a[i] == min_len) {
				hash_dp_init(h_node, hash_pos, start_a, len_a, i, (i-1)*hash_step/*ref offset*/, head, hash_len, hash_step, MIN_FLAG, split_offset, 1/*bound_flag*/, bound_len);
				min_exist = 1;
			}
			else hash_dp_init(h_node, hash_pos, start_a, len_a, i, (i-1)*hash_step/*ref offset*/, head, hash_len, hash_step, MULTI_FLAG, split_offset, 1/*bound_flag*/, bound_len);
		}
	}
	//dp update and backtrack
	{
		//min update and backtrack
		if (min_exist) {
			//min extend
#ifndef __NEW__
			for (i = 1; i <= hash_seed_n; ++i) {
				if (len_a[i] > min_len)
					hash_min_extend(h_node, len_a, i, hash_seed_n+2, min_len, MIN_FLAG, hash_len, hash_step, split_offset);
			}
#endif
#ifdef __NEW__
			if (_head) hash_min_extend(h_node, len_a, head.x, head.y, hash_seed_n+2, min_len, MIN_FLAG, hash_len, hash_step, split_offset);
			for (i = 1; i <= hash_seed_n; ++i) {
				if (len_a[i] <= min_len && len_a[i] > 0) {
					for (j = 0; j < len_a[i]; ++j) 
						hash_min_extend(h_node, len_a, i, j, hash_seed_n+2, min_len, MIN_FLAG, hash_len, hash_step, split_offset);
				}
			}
			if (_tail) hash_min_extend(h_node, len_a, tail.x, tail.y, hash_seed_n+2, min_len, MIN_FLAG, split_offset, hash_len, hash_step);
#endif
			//min update
			for (i = 2; i <= hash_seed_n; ++i) {
				for (j = 0; j < len_a[i]; ++j) {
					if (h_node[i][j].dp_flag == MIN_FLAG)
						hash_dp_update(h_node, len_a, i, j, 1/*update start pos*/, hash_len, hash_step, MIN_FLAG, split_offset, bound_flag, bound_len);
				}
			}
			hash_dp_update(h_node, len_a, tail.x, tail.y, 1/*update start pos*/, hash_len, hash_step, MIN_FLAG, split_offset, bound_flag, bound_len);
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
					if ((right.x == tail.x && _tail == 0 )|| (left.x == head.x && _head == 0)) bound_flag = 0;
					else bound_flag = 1;
					mini_len = mini_hash_main_line(h_node, hash_pos, start_a, len_a, hash_len, hash_step, left, right, _line, split_offset, bound_flag, bound_len);
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
		else	
		{
			for (i = 2; i <= hash_seed_n; ++i) {
				for (j = 0; j < len_a[i]; ++j) {
					if (h_node[i][j].dp_flag == MULTI_FLAG)
						hash_dp_update(h_node, len_a, i, j, 1/*update start pos*/, hash_len, hash_step, MULTI_FLAG, split_offset, bound_flag, bound_len);
				}
			}
			hash_dp_update(h_node, len_a, tail.x, tail.y, 1/*update start pos*/, hash_len, hash_step, MULTI_FLAG, split_offset, bound_flag, bound_len);
			node_i = h_node[tail.x][tail.y].node_n-1;

			int last_x, last_y;
			i = h_node[tail.x][tail.y].from.x;
			j = h_node[tail.x][tail.y].from.y;
			while (i != head.x)
			{
				if (node_i < 0) { fprintf(stderr, "[hash main line] bug: node_i < 0. %d %d\n", _head, _tail); exit(-1); }
				line[node_i].x = i;
				line[node_i].y = j;
				--node_i;
				last_x = h_node[i][j].from.x;
				last_y = h_node[i][j].from.y;
				i = last_x;
				j = last_y;
			}
			if (node_i >= 0) {fprintf(stderr, "[hash main line] bug: node_i >= 0.\n"); exit(-1);}
			return h_node[tail.x][tail.y].node_n;
		}
	}
}


//return: length of node-path
/*int hash_main_line(int *hash_pos, int *start_a, int *len_a, int ref_len, int read_len, int hash_seed_n, int hash_len, int hash_step, hash_dp_node **h_node, line_node *line, int _head, int _tail)
{
	int i, j;
	int node_i;
	int min_len = HASH_MIN_LEN, dp_flag = 0, flag=0;
	line_node head, tail;
	//dp init
		//						from, ref_i, offset, score, node_n, match_flag, dp_flag
		if (_head == 1) 
		{
			head = (line_node){0, 0}; 
			h_node[0][0] = (hash_dp_node){(line_node){-1,0}, 0-hash_len, 0, 0, 0, F_MATCH, MIN_FLAG};
		}
		else 
		{
			head = (line_node){-1, 0};
			h_node[0][0].dp_flag = 0 - MIN_FLAG;
		}
		if (_tail == 1) 
		{
			tail = (line_node){hash_seed_n+1, 0}; 
			h_node[hash_seed_n+1][0] = (hash_dp_node){head, ref_len, read_len-ref_len, 0, 0, F_UNMATCH, MIN_FLAG};
		}
		else 
		{
			tail = (line_node){-1, 0};
			h_node[hash_seed_n+1][0].dp_flag = 0 - MIN_FLAG;
		}

		//for min_len, dp_flag = 1, else 0
		for (i = 1; i <= hash_seed_n; ++i)
		{
			if (len_a[i] == min_len)
			{
				dp_flag = MIN_FLAG;
				flag = 1;
			}
			else dp_flag = MULTI_FLAG;
			hash_dp_init(h_node, hash_pos, start_a, len_a, i, (i-1)*hash_step, head, hash_len, dp_flag);
		}
#ifdef __DEBUG__    
    printf("init:\n");
    for (i = 0; i <= hash_seed_n+1; ++i)
    {
        for (j = 0; j < len_a[i]; ++j)
        {
            printf("%d %d: from %d %d ref_i %d offset %d node_n %d score %d dp_flag %d\n",i,j, h_node[i][j].from.x, h_node[i][j].from.y, h_node[i][j].ref_i, h_node[i][j].offset, h_node[i][j].node_n, h_node[i][j].score, h_node[i][j].dp_flag);
        }
    }
#endif

	if (flag) {	//min_len aln exist
		//min-aln match-extend
#ifndef __NEW__
		for (i = 1; i < hash_seed_n+1; ++i)
		{
			if (len_a[i] > min_len)
				hash_min_extend(h_node, len_a, i, hash_seed_n+2, min_len, MIN_FLAG);
		}
#endif
#ifdef __NEW__
        for (i = 0; i <= hash_seed_n; ++i)
        {
            if (len_a[i] <= min_len)
            {
                for (j = 0; j < len_a[i]; ++j)
                {
                    hash_min_extend(h_node, len_a, i, j, hash_seed_n+2, min_len, MIN_FLAG, hash_len);
                }
            }
        }
		if (_tail == 1)
			hash_min_extend(h_node, len_a, hash_seed_n+1, 0, hash_seed_n+2, min_len, MIN_FLAG, hash_len);
#endif
#ifdef __DEBUG__    
        printf("extend:\n");
        for (i = 0; i <= hash_seed_n+1; ++i)
        {
            for (j = 0; j < len_a[i]; ++j)
            {
                printf("%d %d: from %d %d ref_i %d offset %d node_n %d score %d dp_flag %d\n",i,j, h_node[i][j].from.x, h_node[i][j].from.y, h_node[i][j].ref_i, h_node[i][j].offset, h_node[i][j].node_n, h_node[i][j].score, h_node[i][j].dp_flag);
            }
        }
#endif
        //dp-update for min-aln match-line nodes, whose dp_flag = MIN_FLAG 
		int start = 1;
		for (i = 2; i <= hash_seed_n; ++i)
		{
			for (j = 0; j < len_a[i]; ++j)
			{
				if (h_node[i][j].dp_flag == MIN_FLAG)
					hash_dp_update(h_node, len_a, i, j, start, hash_len, MIN_FLAG);
			}
		}
		if (_tail == 1) hash_dp_update(h_node, len_a, hash_seed_n+1, 0, start, hash_len, MIN_FLAG);
#ifdef __DEBUG__    
        printf("update:\n");
        for (i = 0; i <= hash_seed_n+1; ++i)
        {
            for (j = 0; j < len_a[i]; ++j)
            {
                printf("%d %d: from %d %d ref_i %d offset %d node_n %d score %d dp_flag %d\n",i,j, h_node[i][j].from.x, h_node[i][j].from.y, h_node[i][j].ref_i, h_node[i][j].offset, h_node[i][j].node_n, h_node[i][j].score, h_node[i][j].dp_flag);
            }
        }
#endif

		//dp-update for the blank remaining, backtrack for existing blank
		//line: reverse firstly, convert to be forward
		int mini_len;
		line_node *_line = (line_node*)malloc(hash_seed_n * sizeof(line_node));;
		line_node right, left;
		node_i = 0;

		right = tail;
		left = h_node[tail.x][tail.y].from;
		while (1)
		{
			if (h_node[right.x][right.y].match_flag != F_MATCH && left.x < right.x-1)// there is a seed left
			//if (h_node[left.x][left.y].ref_i 
			{
				mini_len = mini_hash_main_line(h_node, hash_pos, start_a, len_a, hash_len, hash_step, left, right, _line);
				//add mini-line to the main-line, _line is forward, but the line is reverse
				for (i = 0; i < mini_len; ++i)
					line[node_i++] = _line[mini_len - 1 - i];
			}
			if (left.x == head.x)
				break;
			line[node_i++] = left;
			right = left;
			left = h_node[right.x][right.y].from;
		}
		free(_line);
		//convert line to be forward
		line_node tmp;
		for (i = 0; i < node_i/2; ++i)
		{
			tmp = line[i];
			line[i] = line[node_i-1-i];
			line[node_i-1-i] = tmp;	
		}
		return node_i;
	} else { 
		//dp-update for the whole seeds' aln result
        //init again? XXX
		for (i = 2; i < hash_seed_n+1; ++i)
		{
			for (j = 0; j < len_a[i]; ++j)
            {
                if (h_node[i][j].dp_flag == MULTI_FLAG)
                {
                    hash_dp_update(h_node, len_a, i, j, head.x+1, hash_len, MULTI_FLAG);
                }
            }
		}
		hash_dp_update(h_node, len_a, tail.x, tail.y, head.x+1, hash_len, MULTI_FLAG);
#ifdef __DEBUG__
		printf("whole-update\n");
        for (i = 0; i <= hash_seed_n+1; ++i)
        {
            for (j = 0; j < len_a[i]; ++j)
            {
                printf("%d %d: from %d %d ref_i %d offset %d node_n %d score %d dp_flag %d\n",i,j, h_node[i][j].from.x, h_node[i][j].from.y, h_node[i][j].ref_i, h_node[i][j].offset, h_node[i][j].node_n, h_node[i][j].score, h_node[i][j].dp_flag);
            }
        }
#endif

		//line_node right, left;
		node_i = h_node[tail.x][tail.y].node_n-1;
		int last_x, last_y;
		//right = tail;
		i = h_node[tail.x][tail.y].from.x;
		j = h_node[tail.x][tail.y].from.y;
		//left = h_node[right.x][right.y].from;
		while (i != head.x)
		{
			if (node_i < 0) { fprintf(stderr, "[hash main line] bug: node_i < 0.\n"); exit(-1); }
			line[node_i].x = i;
			line[node_i].y = j;
			--node_i;
			last_x = h_node[i][j].from.x;
			last_y = h_node[i][j].from.y;
			i = last_x;
			j = last_y;
		}
		//while (left.x != head.x)
		//{
		//	if (node_i < 0) 
        //    {
        //        fprintf(stderr, "[hash main line] bug: node_i < 0.\n"); 
        //        exit(-1);
        //    }
		//	line[node_i] = left;
		//	node_i--;
		//	right = left;
		//	left = h_node[right.x][right.y].from;
		//}
		if (node_i >= 0) {fprintf(stderr, "[hash main line] bug: node_i >= 0.\n"); exit(-1);}
		return h_node[tail.x][tail.y].node_n;
	}
}*/

/*
int hash_split_map(ucigar32_t **split_cigar, uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, int hash_len, int key_len, ucigar32_t *hash_num, uint64_t **hash_node, cigar32_t *hash_pos)
{
	int i, read_i;//q_key_int, q_kmer_int, 
	uint8_t *ref_query;
	int cur_ref_bound, new_ref_bound, cur_read_bound, new_read_bound;
	int q_key_int, q_kmer_int, node_i;

	//add (ref_len, read_len)
	int *start_a = (int*)malloc((ref_len-hash_len+1+1) * sizeof(int));
	int *len_a = (int*)malloc((ref_len-hash_len+1+1) * sizeof(int));
	if (start_a == NULL || len_a == NULL) {fprintf(stderr, "[hash_split_map] Not enougy memory.\n");}

	for (i = 0; i <= ref_len - hash_len; ++i)
	{
		ref_query = ref_seq + i;
		hash_calcu(&q_key_int, &q_kmer_int, ref_query, hash_len, key_len);
		if (hash_hit(hash_num, hash_node, &node_i, q_key_int, q_kmer_int, hash_len-key_len) == 1)
		{
			start_a[i] = ((hash_node[q_key_int][node_i] & 0xffff0000) >> 16);
			len_a[i] = (hash_node[q_key_int][node_i] & 0xffff);
		}
		else start_a[i] = len_a[i] = 0;	//XXX
	}

	
	   //for (i = 0; i <= ref_len - hash_len; ++i) 
	   //{
	   //int j;
	   //fprintf(stdout, "%d %d -> ", i, len_a[i]);
	   //for (j = 0; j < len_a[i]; ++j)
	   //fprintf(stdout, "%d ", hash_pos[start_a[i]+j]);
	   //fprintf(stdout,"\n");
	   //}
	//main line
	int **line, *line_end, m_i;
	hash_DP_node *h_node;

	line = (int**)malloc(((ref_len-hash_len+1)*(read_len-hash_len+1)+2) * sizeof(int*));
	h_node = (hash_DP_node*)malloc(((ref_len-hash_len+1)*(read_len-hash_len+1)+2) * sizeof(hash_DP_node));
	for (i = 0; i < (ref_len-hash_len+1)*(read_len-hash_len+1)+2; ++i) {
		line[i] = (int*)malloc((ref_len-hash_len+1+1) * sizeof(int));
		//XXX 100 * sizeof
		h_node[i].line.line_i = (int32_t*)malloc(100 * sizeof(int32_t));
	}
	line_end = (int*)malloc(((ref_len-hash_len+1)*(read_len-hash_len+1)+2) * sizeof(int));
	if (line == NULL || h_node == NULL || line_end == NULL)
	{
		fprintf(stderr, "[hash split map]Not enough memory.\n");
		exit(-1);
	}

	//line[m_i][]: [(offset, start, len),(),()...]
	m_i = hash_main_line(hash_pos, start_a, len_a, ref_len, read_len, hash_len, line, line_end, h_node);
	if (line_end[m_i] == 0) {fprintf(stderr, "[hash_split_map] no hash-node is uniquely aligned to the read.\n"); return 0; }

	//fill blank with SV and SW, return the WHOLE cigar
	int _q_len, _t_len, _clen, _b_w, _score; 
	int split_clen=0;
	int line_off;
	uint32_t *_cigar=NULL;
	uint32_t *g_cigar;
	g_cigar = (uint32_t*)malloc(100 * sizeof(uint32_t));

	//1. fix the region between left bound and first line
	if (line[m_i][1] !=0 && (line[m_i][1] + line[m_i][0]) != 0)//blank exists
	{
		//split ref
		_t_len = line[m_i][1]; 
		//split read
		_q_len = line[m_i][1] + line[m_i][0];
		_b_w = abs(_t_len - _q_len);
		_score = ksw_global(_q_len, read_seq, _t_len, ref_seq, 5, &bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
		//fprintf(stdout, "before(%d %d): ", _t_len, _q_len);printcigar(_cigar, _clen);fprintf(stdout, "\n");
		_push_cigar(split_cigar, &split_clen, _cigar, _clen);
	}
	else//add cigar
	{
		//         ref_left, read_left, ref_rigth,    read_right
		make_indel_cigar(-1, -1,        line[m_i][1], line[m_i][1]+line[m_i][0], &_clen, &g_cigar);
		//fprintf(stdout, "before: ");printcigar(g_cigar, _clen);fprintf(stdout, "\n");
		_push_cigar(split_cigar, &split_clen, g_cigar, _clen);
	}
	//2. fix the region betweed each lines
	for (i = 0; i < line_end[m_i]; ++i)
	{
		//generate cigar for every main-lines
		//generate_mem_cigar(,);
		if ((line[m_i][i*3+2]+hash_len-1)==0) _clen = 0;
		else {g_cigar[0] = ((line[m_i][i*3+2]+hash_len-1) << 4) | CMATCH; _clen = 1;}
		//fprintf(stdout, "line: ");printcigar(g_cigar, _clen);fprintf(stdout, "\n");
		_push_cigar(split_cigar, &split_clen, g_cigar, _clen);

		//check blank
		if (i == (line_end[m_i]-1)) break;
		if ((line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2) < (line[m_i][(i+1)*3+1]-1) 
				&& (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]) < (line[m_i][(i+1)*3+1]+line[m_i][(i+1)*3]-1))	//blank exists
		{
			_t_len = (line[m_i][(i+1)*3+1]-1) - (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2);
			_q_len = (line[m_i][(i+1)*3+1]+line[m_i][(i+1)*3]-1) - (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]); 	//blank exists
			_b_w = abs(_t_len - _q_len);
			_score = ksw_global(_q_len, read_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len+line[m_i][i*3])-1, _t_len, ref_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-1), 5, &bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
			//fprintf(stdout, "blank(%d %d): ", _t_len, _q_len);printcigar(_cigar, _clen);fprintf(stdout, "\n");
			_push_cigar(split_cigar, &split_clen, _cigar, _clen);
		}
		else //add cigar
		{
			line_off = make_indel_cigar((line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2), line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3], line[m_i][(i+1)*3+1], line[m_i][(i+1)*3+1]+line[m_i][(i+1)*3], &_clen, &g_cigar);
			//for overlaps between lines
			line[m_i][(i+1)*3+1] -= line_off;
			line[m_i][(i+1)*3+2] += line_off;
			//fprintf(stdout, "sv: ");printcigar(g_cigar, _clen);fprintf(stdout, "\n");
			_push_cigar(split_cigar, &split_clen, g_cigar, _clen);
		}
	}
	//3. fix the region between last line and right bound	//i = line_end[m_i]-1
	if ((line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2) < (ref_len - 1) && (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]) < (read_len-1))	//blank exists
	{
		_t_len = (ref_len-1) - (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2);
		_q_len = (read_len-1) - (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]);
		_b_w = abs(_t_len - _q_len);
		_score = ksw_global(_q_len, read_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len+line[m_i][i*3]-1), _t_len, ref_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-1), 5, &bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
		//_score = ksw_global(_t_len, ref_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-1), _q_len, read_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len+line[m_i][i*3]-1), 5, &bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
		//fprintf(stdout, "last(%d %d): ", _t_len, _q_len);printcigar(_cigar, _clen);fprintf(stdout, "\n");
		_push_cigar(split_cigar, &split_clen, _cigar, _clen);
	}
	else //add cigar
	{
		make_indel_cigar((line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2), (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]), ref_len, read_len, &_clen, &g_cigar);
		//fprintf(stdout, "last: ");printcigar(g_cigar, _clen);fprintf(stdout, "\n");
		_push_cigar(split_cigar, &split_clen, g_cigar, _clen);
	}

	if (_cigar != NULL) free(_cigar); 
	free(g_cigar);

	free(start_a); free(len_a);	
	for (i = 0; i < (ref_len-hash_len+1)*(read_len-hash_len+1)+2; ++i) 
	{ 
		free(line[i]); 
		free(h_node[i].line.line_i);
	}
	free(line); free(line_end);
	free(h_node);
	return split_clen;
}
*/
int hash_split_map(cigar32_t **split_cigar, int *split_clen, int *split_m,
				   uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
				   int hash_len, int hash_step, int key_len, int split_len,
				   uint32_t *hash_num, uint64_t **hash_node, int **hash_node_num, int32_t *hash_pos, 
				   int _head, int _tail)
{
	int i, res=0;   //res: 0x01 -> pull trigger & cut cigar/ 0x10 -> split result
	uint8_t *ref_query;
	int q_key_int, q_kmer_int, node_i;

    //un-overlap XXX
    int hash_seed_n = (ref_len-hash_len)/hash_step + 1;
	int *start_a = (int*)malloc((hash_seed_n + 2) * sizeof(int));	//2: for head and tail
	int *len_a = (int*)malloc((hash_seed_n + 2) * sizeof(int));
	if (start_a == NULL || len_a == NULL) {
        fprintf(stderr, "[hash_split_map] Not enougy memory.(ref_len %d, read_len %d)\n", ref_len, read_len); 
        exit(-1);
    }

	(*split_clen) = 0;

	//hash map, restore hash map result
	{
		len_a[0] = 1;		//for head node
		for (i = 0; i <= ref_len - hash_len; i+=hash_step)
		{
			ref_query = ref_seq + i;
			if (hash_calcu(&q_key_int, &q_kmer_int, ref_query, hash_len, key_len) == 1)	// hash 'N'
			{
				start_a[i/hash_step+1] = 0;
				len_a[i/hash_step+1] = 0;
			} else {
				if (hash_hit(hash_num, hash_node, &node_i, q_key_int, q_kmer_int, hash_len-key_len) == 1)
				{
					start_a[i/hash_step + 1] = ((hash_node[q_key_int][node_i] & 0xffffffff));	//for hit seeds
                    // XXX hash_max_hits
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
		if (line == NULL || h_node == NULL) {fprintf(stderr, "[hash_split_map] Not enougy memory.\n"); exit(-1);} 

		m_len = hash_main_line(hash_pos, start_a, len_a, ref_len, read_len, hash_seed_n, hash_len, hash_step, h_node, line, _head, _tail);
	}

	int _q_len, _t_len, _clen, _b_w, _score; 
	cigar32_t *_cigar=0;
    if (m_len > 0)
	{
		//fill blank with generated SV and SW, return the whole cigar
			cigar32_t *g_cigar;
			g_cigar = (cigar32_t*)malloc(sizeof(cigar32_t));
			int tail_in, head_in;
			//XXX hash read len == 0 XXX
		//1. fix the region between left bound and first line
			if (_head)
			{
				if (h_node[line[0].x][line[0].y].ref_i != 0 && (h_node[line[0].x][line[0].y].ref_i + h_node[line[0].x][line[0].y].offset) != 0)	//blank exists
				{
					tail_in = hash_len/2;
					_t_len = h_node[line[0].x][line[0].y].ref_i+tail_in;
					_q_len = h_node[line[0].x][line[0].y].ref_i + h_node[line[0].x][line[0].y].offset+tail_in;
					_b_w = abs(_t_len - _q_len)+3;
					_score = ksw_global(_q_len, read_seq, _t_len, ref_seq, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
					//XXX plus hash-len all the time, OR just when the score < 0 ???
					/*if (_score < 0)
					{
						int j;
						printf("hash old ref:\t"); for (j = 0; j < _t_len; ++j) printf("%c", "ACGT"[(ref_seq)[j]]); printf("\n");
						printf("hash old read:\t"); for(j = 0; j < _q_len; ++j) printf("%c", "ACGT"[(read_seq)[j]]); printf("\n");
						printf("hash old head cigar:\t");
						printcigar(_cigar, _clen);
						k = 0;
						for (j = 0; j < _clen; ++j)
						{
							if ((_cigar[j] & 0xf) != CMATCH)
								k++;
						}
						fprintf(stdout, "\t%d\t%d\n", k, _clen);

					
						//	_t_len = h_node[line[0].x][line[0].y].ref_i+hash_len;
						//	_q_len = h_node[line[0].x][line[0].y].ref_i + h_node[line[0].x][line[0].y].offset+hash_len;
						//	_b_w = abs(_t_len - _q_len);
						//	_score = ksw_global(_q_len, read_seq, _t_len, ref_seq, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
						//	printf("head hash new ref:\t"); for (j = 0; j < _t_len; ++j) printf("%c", "ACGT"[(ref_seq)[j]]); printf("\n");
						//	printf("hash new read:\t"); for(j = 0; j < _q_len; ++j) printf("%c", "ACGT"[(read_seq)[j]]); printf("\n");
						//	printf("hash new head cigar:\t");
						//	printcigar(_cigar, _clen);
						//	fprintf(stdout, "\n");
					}*/
					//XXX
					if ((_clen>>1) > 10 && _score < 0 && _q_len > 100 && _t_len > 100)
					{
						res |= 1;
                        //free(_cigar);
                        cigar32_t *k_cigar=0;
                        int k_clen = 0;
						if (ksw_both_extend(_q_len, read_seq, _t_len, ref_seq, 5, bwasw_sc_mat, 5, 2, _b_w, hash_len*bwasw_sc_mat[0], hash_len*bwasw_sc_mat[0], &k_clen, &k_cigar) == 1)
                        {
                            free(_cigar);
                            _cigar = k_cigar;
                            _clen = k_clen;
                        }
                    }
					_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
                    free(_cigar);
				}
				else//no blank, add SV cigar
				{	
					//         ref_left, read_left, ref_right,           read_right
					make_indel_cigar(-1, -1,  h_node[line[0].x][line[0].y].ref_i, h_node[line[0].x][line[0].y].ref_i+h_node[line[0].x][line[0].y].offset, &_clen, &g_cigar, split_len, &res);
					_push_cigar(split_cigar, split_clen, split_m, g_cigar, _clen);
					tail_in = 0;
				}
			}
			else tail_in = 0;
		//2. fix the region betweed each match-lines
			int start_i;
			int overlap = 0;	//F_UNMATCH seeds' overlap
			start_i = 0;
			//XXX here only one node
			for (i = 0; i < m_len; ++i)
			{
				if (i == m_len-1 || h_node[line[i+1].x][line[i+1].y].match_flag != F_MATCH)
				{
					//start -> i
					//tail_in
					g_cigar[0] = (h_node[line[i].x][line[i].y].ref_i - h_node[line[start_i].x][line[start_i].y].ref_i + hash_len - tail_in - overlap) << 4 | CMATCH;
					_clen = 1;
					_push_cigar(split_cigar, split_clen, split_m, g_cigar, _clen);
                    //XXX
                    if (g_cigar[0] >>4 > 10000) { fprintf(stderr, "match both tail:%d overlap: %d\t", tail_in, overlap); printcigar(stderr, g_cigar, _clen); printf("\n"); }

					if (i == m_len-1) break;
					//check blank OR generate sv cigar 
					if (h_node[line[i].x][line[i].y].ref_i + hash_len  < h_node[line[i+1].x][line[i+1].y].ref_i 
					 && h_node[line[i].x][line[i].y].ref_i + hash_len + h_node[line[i].x][line[i].y].offset < h_node[line[i+1].x][line[i+1].y].ref_i + h_node[line[i+1].x][line[i+1].y].offset)
					{	//blank exists
						//split_cigar head_in
						if (((*split_cigar)[*split_clen-1] & 0xf) != CMATCH) head_in = 0;
						else 
						{
							if (((*split_cigar)[*split_clen-1] >> 4) > hash_len/2)
							{
								(*split_cigar)[*split_clen-1] -= ((hash_len/2)<<4 | CMATCH);
								//printf("split both "); printcigar(*split_cigar, *split_clen); printf("\n");
								head_in = hash_len/2;
							}
							else
							{
								head_in = ((*split_cigar)[*split_clen-1] >> 4);
								--(*split_clen);
							}
						}
						tail_in = hash_len/2;
						_t_len = h_node[line[i+1].x][line[i+1].y].ref_i - h_node[line[i].x][line[i].y].ref_i - hash_len + head_in + tail_in;
						_q_len = _t_len + h_node[line[i+1].x][line[i+1].y].offset - h_node[line[i].x][line[i].y].offset;
						_b_w = abs(_t_len - _q_len)+3;
						_score = ksw_global(_q_len, read_seq+h_node[line[i].x][line[i].y].ref_i+hash_len+h_node[line[i].x][line[i].y].offset - head_in, _t_len, ref_seq+h_node[line[i].x][line[i].y].ref_i+hash_len - head_in, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
						
						if ((_clen>>1) > 10 && _score < 0 && _q_len > 100 && _t_len > 100) {
							res |= 1;
                            cigar32_t *k_cigar=0; int k_clen=0;
                            if (ksw_both_extend(_q_len, read_seq+h_node[line[i].x][line[i].y].ref_i+hash_len+h_node[line[i].x][line[i].y].offset - head_in, _t_len, ref_seq+h_node[line[i].x][line[i].y].ref_i+hash_len - head_in, 5, bwasw_sc_mat, 5, 2, _b_w, hash_len*bwasw_sc_mat[0], hash_len*bwasw_sc_mat[0], &k_clen, &k_cigar) == 1)
                            {
                                free(_cigar);
                                _cigar = k_cigar;
                                _clen = k_clen;
                            }
                        }
                        _push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
                        overlap = 0;
						free(_cigar);
					}
					else//no blank, add SV cigar
					{
						tail_in = 0;
						overlap = make_indel_cigar(h_node[line[i].x][line[i].y].ref_i+hash_len-1, h_node[line[i].x][line[i].y].ref_i+hash_len-1+h_node[line[i].x][line[i].y].offset, h_node[line[i+1].x][line[i+1].y].ref_i, h_node[line[i+1].x][line[i+1].y].ref_i+h_node[line[i+1].x][line[i+1].y].offset, &_clen, &g_cigar, split_len, &res);
						_push_cigar(split_cigar, split_clen, split_m, g_cigar, _clen);
					}
					start_i = i+1;
				}
				else continue;
			}
		//3. fix the region between last line and right bound	//i = line_end[m_i]-1
			if (_tail)
			{
				//XXX end node, hash_step
				if (h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len < ref_len && h_node[line[m_len-1].x][line[m_len-1].y].ref_i+h_node[line[m_len-1].x][line[m_len-1].y].offset + hash_len < read_len)//blank exists
				{
					if (((*split_cigar)[*split_clen-1] & 0xf) != CMATCH) head_in = 0;
					else 
					{
						if (((*split_cigar)[*split_clen-1] >> 4) > hash_len/2)
						{
							(*split_cigar)[*split_clen-1] -= ((hash_len/2)<<4 | CMATCH);
							//printf("tail both "); printcigar(*split_cigar, *split_clen); printf("\n");
							head_in = hash_len/2;
						}
						else
						{
							head_in = (*split_cigar)[*split_clen-1] >> 4;
							--(*split_clen);
						}
					}
					_t_len = ref_len - h_node[line[m_len-1].x][line[m_len-1].y].ref_i - hash_len + head_in;
					_q_len = read_len - h_node[line[m_len-1].x][line[m_len-1].y].ref_i - hash_len - h_node[line[m_len-1].x][line[m_len-1].y].offset + head_in;
					_b_w = abs(_t_len - _q_len)+3;
                    _score = ksw_global(_q_len, read_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+h_node[line[m_len-1].x][line[m_len-1].y].offset+hash_len - head_in, _t_len, ref_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len - head_in, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);

                    if ((_clen>>1) > 10 && _score < 0 && _q_len > 100 && _t_len > 100)
                    {
                        res |= 1;
                        cigar32_t *k_cigar=0;
                        int k_clen = 0;
                        if (ksw_both_extend(_q_len, read_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+h_node[line[m_len-1].x][line[m_len-1].y].offset+hash_len - head_in, _t_len, ref_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len - head_in, 5, bwasw_sc_mat, 5, 2, _b_w, hash_len*bwasw_sc_mat[0], hash_len*bwasw_sc_mat[0], &k_clen, &k_cigar) == 1)
                        {
                            free(_cigar);
                            _cigar = k_cigar;
                            _clen = k_clen;
                        }
                    }
					_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
					free(_cigar);
				}
				else //no blank, add SV cigar
				{
					make_indel_cigar(h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len-1, h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len-1+h_node[line[m_len-1].x][line[m_len-1].y].offset, ref_len, read_len, &_clen, &g_cigar, split_len, &res);
					_push_cigar(split_cigar, split_clen, split_m, g_cigar, _clen);
				}
			}
		//free variables
		free(g_cigar);
	}
	else	//no hash-dp line nodes exist
	{
		if (_head && _tail)
		{
			_t_len = ref_len; _q_len = read_len;
			_b_w = abs(_t_len - _q_len)+3;
			_score = ksw_global(_q_len, read_seq, _t_len, ref_seq, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
			if ((_clen>>1) > 10 && _score < 0 && _q_len > 100 && _t_len > 100)
			{
				res |= 1;
                cigar32_t *k_cigar=0;
                int k_clen = 0;
                if (ksw_both_extend(_q_len, read_seq, _t_len, ref_seq, 5, bwasw_sc_mat, 5, 2, _b_w, hash_len*bwasw_sc_mat[0], hash_len*bwasw_sc_mat[0], &k_clen, &k_cigar) == 1)
                {
                    free(_cigar);
                    _cigar = k_cigar;
                    _clen = k_clen;
                }
            }
			_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
			free(_cigar);
		}
        //_head or _tail -> left or right
	}
    //_FREE_V:
	for (i = 0; i < hash_seed_n+2; ++i) free(h_node[i]);
	free(h_node); free(line); 
	free(start_a); free(len_a);
	return res;
}

//XXX hash_node: use struct?
//insertion length = read_len - ref_len
// return: 0x00->normal/ 0x01->trigger/ 0x10->split result
int split_indel_map(cigar32_t **res_cigar, int *res_len, int *res_m,
					 uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, 
					 int64_t ref_offset, 
					 int hash_len, int hash_step, int split_len,
					 uint32_t **hash_num, uint64_t ***hash_node, 
					 int key_len, int hash_size)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t)); int **hash_node_num;
	//init of hash table
	init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, key_len, hash_size);

	int res;
	//for multi-breakpoints
	res = hash_split_map(res_cigar, res_len, res_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, key_len, split_len, *hash_num, *hash_node, hash_node_num, hash_pos, 1, 1);
	free(hash_pos);
    int i; for (i = 0; i < hash_size; ++i) free(hash_node_num[i]);
    free(hash_node_num);

	return res;
}

int hash_right_bound_map(cigar32_t **cigar, int *cigar_len, int *cigar_m,
						 uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
						 uint32_t **hash_num, uint64_t ***hash_node, 
						 int hash_len, int hash_key, int hash_step,
                         int split_len)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	int hash_size = pow(NT_N, hash_key);
    int **hash_node_num;
	if (init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, hash_key, hash_size) != 0)
	{
		(*cigar_len) = 0;
		free(hash_pos);
		return 1;
	}
	int res = 0;
	res |= hash_split_map(cigar, cigar_len, cigar_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, hash_key, split_len, *hash_num, *hash_node, hash_node_num, hash_pos, 1, 0);
	int i, readINcigar=0, refINcigar=0;
	for (i = 0; i < (*cigar_len); ++i)
	{
		if (((*cigar)[i] & 0xf) == CMATCH) readINcigar += ((*cigar)[i] >> 4), refINcigar += ((*cigar)[i] >> 4);
		else if (((*cigar)[i] & 0xf) == CINS || ((*cigar)[i] & 0xf) == CSOFT_CLIP) readINcigar += ((*cigar)[i] >> 4);
		else if (((*cigar)[i] & 0xf) == CDEL || ((*cigar)[i] & 0xf) == CHARD_CLIP) refINcigar += ((*cigar)[i] >> 4);
		else 
        {
            fprintf(stderr, "[hash_right_bound_map] cigar error: ");
            printcigar(stderr, *cigar, *cigar_len);
            fprintf(stderr, "\n");
            exit(-1);
        }
	}
	if (readINcigar < read_len)     
	{
		int qlen = read_len-readINcigar;
		int tlen = read_len -readINcigar + ((hash_len+read_len-readINcigar)*bwasw_sc_mat[0]-5-1)/2;
		tlen = tlen < (ref_len - refINcigar) ? tlen : (ref_len - refINcigar);
		int read_end, ref_end, n_cigar_, m_cigar_;
		cigar32_t *cigar_= NULL;
		//
		//printf("ref:\t");for (i = 0; i < tlen; ++i) printf("%c", "ACGT"[(ref_seq+refINcigar)[i]]);printf("\n");
		//printf("read:\t");for (i =0; i < qlen; ++i) printf("%c", "ACGT"[(read_seq+readINcigar)[i]]);printf("\n");
		res |= ksw_extend_c(qlen, read_seq+readINcigar, tlen , ref_seq+refINcigar, 5, bwasw_sc_mat, 5, 2, abs(read_len-readINcigar-ref_len + refINcigar)+3, hash_len*bwasw_sc_mat[0], 4, &read_end, &ref_end, &n_cigar_, &cigar_, &m_cigar_);
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
						uint32_t **hash_num, uint64_t ***hash_node, 
						int hash_len, int hash_key, int hash_step,
                        int split_len)
{
	int hash_cigar_len=0, hash_cigar_m=100;
	cigar32_t *hash_cigar = (cigar32_t*)malloc(hash_cigar_m * sizeof(cigar32_t));;
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	(*cigar_len) = 0;
	int hash_size = (int)pow((double)NT_N, (double)hash_key);	
    int **hash_node_num;
	if (init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, hash_key, hash_size) != 0) {
		(*cigar_len) = 0;
		free(hash_pos); free(hash_cigar);
		return 1;
	}
	int res = 0;
    res |= hash_split_map(&hash_cigar, &hash_cigar_len, &hash_cigar_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, hash_key, split_len, *hash_num, *hash_node, hash_node_num, hash_pos, 0, 1);
	int i, readINcigar=0, refINcigar=0;
	for (i = 0; i < hash_cigar_len; ++i)
	{
		if ((hash_cigar[i] & 0xf) == CMATCH) readINcigar+=(hash_cigar[i] >> 4), refINcigar+=(hash_cigar[i] >> 4);
		else if ((hash_cigar[i] & 0xf) == CINS || (hash_cigar[i] & 0xf) == CSOFT_CLIP) readINcigar += (hash_cigar[i] >> 4);
		else if ((hash_cigar[i] & 0xf) == CDEL || (hash_cigar[i] & 0xf) == CHARD_CLIP) refINcigar += (hash_cigar[i] >> 4);
		else fprintf(stderr, "[hash_left_bound_map] cigar error: "), printcigar(stderr, hash_cigar, hash_cigar_len), fprintf(stderr, "\n"), exit(-1);
	}
	if (readINcigar < read_len)     //'S' exists
	{
		int qlen = read_len - readINcigar;
		int tlen = read_len - readINcigar + ((hash_len+read_len-readINcigar)*bwasw_sc_mat[0]-5-1)/2;
		tlen = tlen < (ref_len -refINcigar) ? tlen : (ref_len - refINcigar);
		uint8_t *qseq = (uint8_t*)malloc(qlen * sizeof(uint8_t));
		uint8_t *tseq = (uint8_t*)malloc(tlen * sizeof(uint8_t));
		int read_end, ref_end, n_cigar_, m_cigar_;
		cigar32_t *cigar_= 0;
		for (i = 0; i < qlen; ++i) qseq[i] = read_seq[qlen-i-1];
		for (i = 0; i < tlen; ++i) tseq[i] = ref_seq[ref_len-refINcigar-i-1];

		res |= ksw_extend_c(qlen, qseq, tlen, tseq, 5, bwasw_sc_mat, 5, 2, abs(qlen-tlen)+3, hash_len*bwasw_sc_mat[0], 4, &read_end, &ref_end, &n_cigar_, &cigar_, &m_cigar_);

		if (cigar_ != NULL) {
			if (read_end < read_len-readINcigar) _push_cigar1(&cigar_, &n_cigar_, &m_cigar_, (((read_len-readINcigar-read_end) << 4) | CSOFT_CLIP));//'S' exsit
			//invert cigar
			int32_t tmp;
			for (i = 0; i < n_cigar_/2; ++i) tmp = cigar_[i], cigar_[i] = cigar_[n_cigar_-i-1], cigar_[n_cigar_-i-1] = tmp;
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
						uint32_t **hash_num, uint64_t ***hash_node, 
						int hash_len, int hash_key, int hash_step,
                        int split_len)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	int hash_size = (int)pow((double)NT_N, (double)hash_key);	
    int **hash_node_num;
	if (init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_node_num, &hash_pos, hash_key, hash_size) != 0) {
		(*cigar_len = 0);
		free(hash_pos);
		return 1;
	}
	int res = 0;
	res |= hash_split_map(cigar, cigar_len, cigar_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, hash_key, split_len, *hash_num, *hash_node, hash_node_num, hash_pos, 1, 1);
	/*int i, readINcigar=0;
	for (i = 0; i < (*cigar_len); ++i)
	{
		if (((*cigar)[i] & 0xf) == CMATCH || ((*cigar)[i] & 0xf) == CINS)
			readINcigar+=((*cigar)[i] >> 4);
	}
	if (readINcigar < read_len)     //'S' exists
	{
		//LV?
		for (i = (*cigar_len); i > 0; --i)
		{
			(*cigar)[i] = (*cigar)[i-1];
		}
		(*cigar)[0] = (((read_len - readINcigar) << 4) | CSOFT_CLIP);
		(*cigar_len)++;
	}*/
	
	free(hash_pos);
    int i; for (i = 0; i < hash_size; ++i) free(hash_node_num[i]);
    free(hash_node_num);
	return res;
}
