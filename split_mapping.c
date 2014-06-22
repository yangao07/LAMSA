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
	for (i = 0; i < key_len; ++i) 
	{
		if (seed[i] >= (1 << key_len))
			return 1;
		(*key_int) = (*key_int) << 2 | (int)seed[i];
	}
	for (; i < hash_len; ++i) 
	{
		if (seed[i] >= (1 << key_len))
			return 1;
		(*kmer_int) = (*kmer_int) << 2 | (int)seed[i];	
	}
	//if (*key_int > 16) {fprintf(stderr, "[hash init] Error: key_int %d\thash_len: %d\tkey_len: %d\nseed:\t", *key_int, hash_len, key_len); for (i = 0; i < hash_len; ++i) fprintf(stderr, "%d",seed[i]); fprintf(stderr, "\n");exit(-1); }
	
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
				   uint32_t *hash_num, uint64_t ***hash_node, int32_t **hash_pos)
{
	int key_int, kmer_int, node_i;
	
	hash_calcu(&key_int, &kmer_int, hash_seq, hash_len, key_len);
	//node_i = hash_search(&hash_num, hash_node, key_int, kmer_int, hash_len-key_len, &EQUAL_FLAG);
	if (hash_hit(hash_num, *hash_node, &node_i, key_int, kmer_int, hash_len-key_len) != 1)
		fprintf(stderr, "2rd search wrong.\n");
	//printf("i: %d start: %d alr: %d ", seed_i, (((*hash_node)[key_int][hash_i] & 0xff00) >> 8), ((*hash_node)[key_int][hash_i] & 0xff));
	(*hash_pos)[(((*hash_node)[key_int][node_i] & 0xffff0000) >> 16) + ((*hash_node)[key_int][node_i] & 0xffff)] = hash_offset;
	++(*hash_node)[key_int][node_i];
	return 0;
}

	//@function:
	//hash_node:  ----32----||--------16--------|---------16----------       
	//             kmer_int           0               pos_num
	//                                ->
	//hash_node:  ----32----||--------16--------|---------16----------       
	//             kmer_int    pos_start_index   pos_already_in_num(0)
int hash_calcu_pos_start(int hash_size, uint64_t ***hash_node, uint32_t *hash_num)
{
	int i, j, tmp, pos_start=0;
	for (i = 0; i < hash_size; ++i)
	{
		for (j = 0; j < hash_num[i]; ++j)
		{
			//printf("pos_start: %d\n", pos_start);
			tmp = pos_start;
			pos_start += ((*hash_node)[i][j]&0xffff);
			(*hash_node)[i][j] = ((*hash_node)[i][j]&0xffffffff00000000) | (tmp<<16);
		}
	}
	//printf("pos_start: %d\n", pos_start);
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
			  uint32_t **hash_num, uint64_t ***hash_node, int32_t **hash_pos, 
			  int key_len, int hash_size)
{
	int i;
	uint8_t *hash_seq;

	//check 'N', repalce by 'G'
	for (i = 0; i < read_len; ++i) if (read_seq[i] > 3) read_seq[i] = 2;
	for (i = 0; i < hash_size; ++i)
		(*hash_num)[i] = 0;
	for (i = 0; i <= read_len-hash_len; ++i)	//create hash index
	{
		//strncpy(kmer, read_seq,	hash_len);
        hash_seq = read_seq+i;
		init_hash_pos_num(hash_seq, hash_len, key_len, hash_num, hash_node);
	}
	hash_calcu_pos_start(hash_size, hash_node, *hash_num);
	for (i = 0; i <= read_len-hash_len; ++i)
	{
		hash_seq = read_seq+i;
		init_hash_core(hash_seq, i, hash_len, key_len, hash_size, *hash_num, hash_node, hash_pos);
	}
	//printf("hash_pos:\n");
	//for (i = 0; i < read_len; ++i) printf("%d ", (*hash_pos)[i]);
	//printf("\n");
	return 0;
}

/*int hash_exact_map(uint64_t **hash_node, int key_int, int node_i, int32_t *hash_pos, int offset)
{
	int start, len, res, i;

	start = ((hash_node[key_int][node_i] & 0xffff0000) >> 16);
	len = (hash_node[key_int][node_i] & 0xffff);
	if (len < 1) fprintf(stderr, "[hash_map] ERROR: hit-len < 1.\n");
	res = hash_pos[start];
	for (i = 1; i < len; ++i)
	{
		if (abs(hash_pos[start+i] - offset) < abs(res - offset))
			res = hash_pos[start+i];
	}
	return res;
}*/

/*int check_hash(int h_offset, int h_i)
{
    int dis = h_offset - h_i;
    if (dis <= 10 && dis >= -10)
        return F_MATCH;
    else if (dis > 10)
        return F_DELETE;
    else return F_INSERT;
}*/

/*int hash_map(uint32_t *hash_num, uint64_t **hash_node, int32_t *hash_pos, int *r_i, uint8_t *query, int offset, int hash_len, int key_len)
{
	int q_key_int, q_kmer_int, node_i;

	hash_calcu(&q_key_int, &q_kmer_int, query, hash_len, key_len);
    //printf("%d %d ", q_key_int, q_kmer_int);
	//hash-hit is different with hash-search?
	if (hash_hit(hash_num, hash_node, &node_i, q_key_int, q_kmer_int, hash_len-key_len) == 1)//hit
	{
		*r_i = hash_exact_map(hash_node, q_key_int, node_i, hash_pos, offset);
		//printf("ref %d read %d \n", offset, *r_i);
        if (check_hash(offset, *r_i) == F_MATCH)
            return 1;
        else return 2;
	}
	else return 0;
}*/

/*int hash_get_dis(int *hash_pos, int *start_a, int a, int a_i, int b, int b_i, int hash_len, int *con_flag)
{
	if (a == b && a_i == b_i) {*con_flag = F_MATCH; return 0;}
	int exp, act, dis;
	//if (abs(a-b) < 10)
	
	//abs(a-b) >= 10:
	exp = hash_pos[start_a[a]+a_i] + b - a;
	act = hash_pos[start_a[b]+b_i];
	dis = (abs(b-a)/(b-a)) * (act - exp);


	if (dis >= 5 && dis < DEL_THD) *con_flag = F_INSERT;
    else if (dis <= -5 && dis >= (0-(abs(a-b)-hash_len))) *con_flag = F_DELETE;
    else if (dis < 5 && dis > -5) *con_flag = F_MATCH; 
	else *con_flag = F_UNCONNECT;
	//printf("%d", dis);
	return abs(dis);
}*/

//1-mismatch betweed two adjacent nodes is allowed
//@para: bound_flag: 1 -> _head,_tail both are set as 1, 0 -> only one.
//对于每一个节点，依据于已经存在的 split_offset，设定其对应的 split_offset	XXX
int hash_main_dis(int a_i, int a_offset, int b_i, int b_offset, 
				  int hash_len, int hash_step, 
				  int *con_flag, 
				  int split_offset, 	//XXX
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
		/*else
		{
			if (abs(b_i - a_i) <= hash_len + 1)
				*con_flag = F_MATCH;
			else if (abs(b_i - a_i) <= hash_len + hash_len)	//the second 'hash_len' represents the last hash_len bps that can't construct a seed
				*con_flag = F_MISMATCH;
			else *con_flag = F_UNCONNECT;
		}*/
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

/*int hash_new_line(int offset, int hash_i, int **line, int *line_end, int path_n)
{
	//line[i][0]: the offset (read_pos - ref_pos)
	line[path_n][0] = offset;
	line[path_n][1] = hash_i;
	line_end[path_n] = 1;
	return 0;
}*/

/*int hash_copy_line(int **line, int *line_end, int from, int m_i)
{
	//for (i = 1; i <= line_end[from]; ++i)
	//{
	//	line[to][line_end[to]+i-1] = line[from][i];
	//}
	//line_end[to] += line_end[from];

	int i,j;
	for (i = 0; i < line_end[m_i]; ++i)
	{
		if (line[m_i][i*3+1] > line[from][1])
			break;
	}
	for (j = line_end[m_i]; j > i; --j)
	{
		line[m_i][j * 3] = line[m_i][(j-1) * 3];	//offset
		line[m_i][j * 3 + 1] = line[m_i][(j-1) * 3 + 1];	// start
		line[m_i][j * 3 + 2] = line[m_i][(j-1) * 3 + 2];	//num of covered nodes
	}
	line[m_i][i * 3] = line[from][0];	//offset
	line[m_i][i * 3 + 1] = line[from][1];	//start
	line[m_i][i * 3 + 2] = line[from][line_end[from]] - line[from][1] + 1;	//num of covered nodes
	++line_end[m_i];
	return 0;
}*/

/*int hash_add_line(int **line, int *line_end, int line_i, int seed_i)
{
	int i, j;
	for (i = 1; i <= line_end[line_i]; ++i)
	{
		if (line[line_i][i] > seed_i)
		{
			break;
		}
		else if (line[line_i][i] == seed_i)
		{
			fprintf(stderr, "[hash_add_line] error.\n");
			exit(-1);
		}
	}
	for (j = line_end[line_i]+1; j > i; --j)
	{
		line[line_i][j] = line[line_i][j-1];
	}
	line[line_i][j] = seed_i;
	++line_end[line_i];
	
	return 0;
}*/

/*hash_line_t *hash_init_line(int seq_len)
{
	hash_line_t *line;
	line = (hash_line_t*)malloc(sizeof(hash_line_t));
	line->n_blanks = 1;	//XXX initial value=1, MAX value=10.
	line->blank = (hash_blank_t*)malloc(10 * sizeof(hash_blank_t));
	line->blank[0].offset = 0;
	line->blank[0].len = seq_len;

	return line;
}*/

	/*int hash_modify_blank(hash_line_t *h_line, int blank_num, int b_start, int b_end)
	{
		if (blank_num >= 10)	//XXX 10
		{
			fprintf(stderr, "[hash_modify_blank] ERROR: blank number is over the range.\n"); exit(-1);
		}
		h_line->blank[blank_num].offset = b_start;
		h_line->blank[blank_num].len = b_end-b_start+1;
		return 0;
	}*/
	/*int hash_delete_blank(hash_line_t *h_line, int blank_num)
	{
		if (blank_num >= 10)	//XXX 10
		{
			fprintf(stderr, "[hash_delete_blank] ERROR: blank number is over the range.\n"); exit(-1);
		}
		int i;
		for (i = blank_num; i < 9; ++i)	//XXX 9
		{
			h_line->blank[i].offset = h_line->blank[i+1].offset;
			h_line->blank[i].len = h_line->blank[i+1].len;
		}
		--(h_line->n_blanks);
		return 0;
	}*/
	/*int hash_insert_blank(hash_line_t *h_line, int blank_num, int b_start, int b_end)
	{
		if (h_line->n_blanks == 10 || blank_num >= 10)	//XXX 10
		{
			fprintf(stderr, "[hash_insert_blank] ERROR: blank number is over the range.\n"); exit(-1);
		}
		int i;
		for (i = h_line->n_blanks; i > blank_num+1; --i)
		{
			h_line->blank[i].offset = h_line->blank[i-1].offset;
			h_line->blank[i].len = h_line->blank[i-1].len;
		}
		h_line->blank[blank_num+1].offset = b_start;
		h_line->blank[blank_num+1].len = b_end-b_start+1;
		++(h_line->n_blanks);
		return 0;
	}*/
					/*ref_line | read_line*/ /*ref_len, read_len*/
	/*int hash_fill_blank(hash_line_t *h_line, int seq_len, int blank_num, int line_start, int line_len)
	{
		if (blank_num < 0) return 0;
		int line_end = line_start+line_len-1;
		int b_start = h_line->blank[blank_num].offset;
		int b_end = b_start + h_line->blank[blank_num].len - 1;
		
		//insert blank(modify AND insert); delete blank; modify blank;
		if (line_end < b_start || line_start > b_end)
		{
			fprintf(stderr, "[hash_fill_blank] Line do NOT match blank.\n"); 
			fprintf(stderr, "                  blank: %d -> %d\n", b_start, b_end);
			fprintf(stderr, "                  line:  %d %d %d\n", line_start, line_len, line_end);
			exit(-1);
		}
		if (line_start <= b_start)
		{
			if (line_end >= b_end)	//delete blank
			{
				hash_delete_blank(h_line, blank_num);
			}
			else //line_end < b_end: modify blank-> (line_end+1, b_end)
			{
				hash_modify_blank(h_line, blank_num, line_end+1, b_end);
			}
		}
		else	//line_start > b_start
		{
			if (line_end < b_end)	//modify AND insert
			{
				hash_modify_blank(h_line, blank_num, b_start, line_start-1);
				hash_insert_blank(h_line, blank_num, line_end+1, b_end);
			}
			else //line_end >= b_end: modify
			{
				hash_modify_blank(h_line, blank_num, b_start, line_start-1);
			}
		}


		return 0;
	}*/

	//blank_cover: return the covered bps, set blank_i as the index of blank
	//if blank_i equal -1, that indicates that this node-line doesn't cover any blank
	/*int blank_cover(hash_line_t *h_line, int *blank_i, int pos_start, int len)
	{
		int covered_len = 0;
		int i;
		int pos_end = pos_start + len - 1;

		(*blank_i = -1);
		for (i = 0; i < h_line->n_blanks; ++i)
		{
			if (h_line->blank[i].offset > pos_end || (h_line->blank[i].offset + h_line->blank[i].len - 1) < pos_start)
				continue;
			else	// have a covered part
			{
				covered_len = (pos_end < (h_line->blank[i].offset + h_line->blank[i].len-1) ? pos_end : (h_line->blank[i].offset + h_line->blank[i].len-1)) - (pos_start > h_line->blank[i].offset ? pos_start : h_line->blank[i].offset) + 1;
				(*blank_i) = i;
				break;
			}
		}
		return covered_len;
	}*/

/*int hash_line_bound(int **line, int *line_end, 
					int *start_a, int *hash_pos, 
					int line_i, int hash_len, 
					int *ref_start, int *ref_end, int *read_start, int *read_end)
{
	int offset = line[line_i][0];

	*ref_start = line[line_i][1];
	*ref_end = line[line_i][line_end[line_i]] + hash_len - 1;
	*read_start = (*ref_start) + offset;
	*read_end = (*ref_end) + offset;

	return 0;
}*/

//pick_the_line: find out the node-line which cover the most bps in blank, set the ref_blank_i and read_blank_i, return the line index.
//if the return value equals -1, it indicates that NO node-line covers any blank.
/*int pick_the_line(int **line, int *line_end, int *copy_flag, int *start_a, int *hash_pos, int hash_len, int path_n, hash_line_t *ref_line, hash_line_t *read_line, int *ref_blank_i, int *read_blank_i)
{
	int i, max_i, max_cover, cover, tmp_ref_blank_i, tmp_read_blank_i;
	int ref_start, ref_end, read_start, read_end;

	max_i = -1;
	max_cover = 0;	//XXX
	for (i = 0; i < path_n; ++i)
	{
		if (copy_flag[i] != 1)
		{
			hash_line_bound(line, line_end, start_a, hash_pos, i, hash_len, &ref_start, &ref_end, &read_start, &read_end);
			//tmp_blank_i == -1  means that the line can't fill any blanks.
			cover = blank_cover(ref_line, &tmp_ref_blank_i, ref_start, ref_end-ref_start+1) //+ blank_cover(read_line, &tmp_read_blank_i, read_start, read_end-read_start+1);
				  * blank_cover(read_line, &tmp_read_blank_i, read_start, read_end-read_start+1)
				  - abs(line[i][0]);	//abs(offset)

			if (cover > max_cover)
			{
				max_i = i;
				*ref_blank_i = tmp_ref_blank_i;
				*read_blank_i = tmp_read_blank_i;
				max_cover = cover;
			}
		}
	}
	return max_i;
}*/

//cut long node-line which has big gap
/*int hash_move_line(int **line, int *line_end, int from, int start, int end, int to)
{
	if (from == to) return 0;

	int i;
	line[to][0] = line[from][0];
	for (i = start; i <= end; ++i)
	{
		line[to][i-start+1] = line[from][i];
	}
	line_end[to] = end - start + 1;
	return 0;
}*/

//XXX add mis-match negative points
/*int hash_cut_line(int **line, int *line_end, int *line_n, int hash_len)
{
	int i, j, tmp_n;

	tmp_n = *line_n;
	for (i = 0; i < tmp_n; ++i)
	{
		if (line_end[i] < 2) continue;
		for (j = line_end[i]; j > 1; --j)
		{
			if ((line[i][j] - line[i][j-1]) > (hash_len+1))
			{
				//move and modify line_end
				hash_move_line(line, line_end, i, j, line_end[i], *line_n);
				line_end[i] = j - 1;
				++(*line_n);
			}
		}
	}
	return 0;
}*/

/*int hash_DP_init(hash_DP_node *h_node, int **line, int *line_end, int l_i, int hash_len)
{
	h_node[l_i].node_i = l_i;
	h_node[l_i].from_i = -1;

	h_node[l_i].line.n_line = 1;
	h_node[l_i].line.line_i[0] = l_i;
	h_node[l_i].line.ref_cover_len = h_node[l_i].line.read_cover_len = line[l_i][line_end[l_i]] - line[l_i][1] + hash_len;
	
		//h_node[l_i].ref_blank = ;
		//h_node[l_i].read_blank = ;
		//h_node[l_i].mis_match;
	return 0;
}*/

/*int hash_update_line(hash_DP_node *h_node, int **line, int *line_end, int node_i, int hash_len)
{
	int i, j, l_i;
	int tmp_ref_cover_plus, tmp_read_cover_plus, max_ref_cover_plus, max_read_cover_plus, max_from_i = -1, max=0;
	int overlap;

	for (i = 0; i < node_i; ++i)
	{
		if (line_end[i] == 0) continue;
		//1. check if it(i) is allowed to be combined with node_i	
		j = h_node[i].line.n_line - 1;	//last line
		l_i = h_node[i].line.line_i[j];
		overlap = (line[node_i][1] > (line[l_i][line_end[l_i]]+hash_len-1)) ? 
			0 : ((line[l_i][line_end[l_i]]+hash_len-1) - line[node_i][1] + 1);
		if (overlap >= (line[node_i][line_end[node_i]] - line[node_i][1] + hash_len))
			continue;
		tmp_ref_cover_plus = h_node[i].line.ref_cover_len - overlap; 
		overlap = ((line[node_i][1]+line[node_i][0]) > (line[l_i][line_end[l_i]]+hash_len-1+line[l_i][0])) ? 
			0 : ((line[l_i][line_end[l_i]]+hash_len-1+line[l_i][0]) - (line[node_i][1]+line[node_i][0]) + 1);
		if (overlap >= (line[node_i][line_end[node_i]] - line[node_i][1] + hash_len))
			continue;
		tmp_read_cover_plus = h_node[i].line.read_cover_len - overlap; 
		//2. check if it(i) is the longest one
		if ((tmp_ref_cover_plus + tmp_read_cover_plus) > max)
		{
			//XXX 
			if ((tmp_ref_cover_plus + tmp_read_cover_plus - max) <= hash_len && max_from_i != -1 && h_node[i].line.n_line > h_node[max_from_i].line.n_line)
				continue;
			max = tmp_ref_cover_plus + tmp_read_cover_plus; 
			max_ref_cover_plus = tmp_ref_cover_plus;
			max_read_cover_plus = tmp_read_cover_plus;
			max_from_i = i;
		}
	}
	if (max_from_i != -1)//update node
	{
		h_node[node_i].from_i = max_from_i;
		h_node[node_i].line.n_line = h_node[max_from_i].line.n_line + 1;
		for (i = 0; i < h_node[max_from_i].line.n_line; ++i)
			h_node[node_i].line.line_i[i] = h_node[max_from_i].line.line_i[i];
		h_node[node_i].line.line_i[i] = node_i;

		h_node[node_i].line.ref_cover_len += max_ref_cover_plus;
		h_node[node_i].line.read_cover_len += max_read_cover_plus;
	}

	return 0;
}*/

/*int hash_swap_line(int **line, int *line_end, int ii, int jj)
{
	int i, tmp;
	if (line_end[ii] > line_end[jj])
	{
		for (i = 0; i <= line_end[jj]; ++i)
		{
			tmp = line[ii][i];
			line[ii][i] = line[jj][i];
			line[jj][i] = tmp;
		}
		for (i = line_end[jj]+1; i <= line_end[ii]; ++i)
		{
			line[jj][i] = line[ii][i];
		}
	}
	else
	{
		for (i = 0; i <= line_end[ii]; ++i)
		{
			tmp = line[jj][i];
			line[jj][i] = line[ii][i];
			line[ii][i] = tmp;
		}
		i = line_end[ii]+1; 
		while(i <= line_end[jj])
		{
			line[ii][i] = line[jj][i];
			++i;
		}
	}
	tmp = line_end[jj];
	line_end[jj] = line_end[ii];
	line_end[ii] = tmp;
	return 0;
}*/


/*int mini_main_line(int *hash_pos, int *start_a, int *len_a, int ref_start, int ref_end, int read_start, int read_end, int hash_len, int **line, int *line_end, hash_DP_node *h_node)
{
    if (((ref_end - ref_start) < hash_len) || ((read_end - read_start) < hash_len))
        return -1;

    int i, j, k, path_n;
    int con_flag, flag;
    int copy_flag[500];

    //generate all node-line
    path_n = 0;
    for (i = ref_start; i <= ref_end-hash_len; ++i)
    {
		for (j = 0; j < len_a[i]; ++j)
		{
			if (hash_pos[start_a[i]+j] <= read_end && (hash_pos[start_a[i]+j]+hash_len-1) >= read_start)
			{
				flag = 0;
				for (k = path_n-1; k >= 0; --k)
				{
					hash_main_dis(, &con_flag);
					if (con_flag == F_MATCH)
					{
						flag = 1;
						line[k][line_end[k]+1] = i;
						++line_end[k];
						break;
					}
				}
				if (flag == 0)
				{
					hash_new_line(hash_pos[start_a[i]+j]-i, i, line, line_end, path_n);
					++path_n;
				}
			}
		}
	}

	if (path_n == 0) return -1;
	hash_cut_line(line, line_end, &path_n, hash_len);
	for (i = 0; i < path_n; ++i) { copy_flag[i] = 0; }
	for (i = 0; i < path_n-1; ++i)
	{
		for (j = i+1; j < path_n; ++j)
		{
			if ((line[j][1] < line[i][1]) || ((line[j][1] == line[i][1]) && (line[j][0] < line[i][0])))
				hash_swap_line(line, line_end, i, j);
		}
	}
	//fprintf(stdout, "all mini line:\nstart:\t");
	for (i = 0; i < path_n; ++i)
	{
		if (line_end[i] != 0)
			fprintf(stdout, "%d ", line[i][1]);
	}fprintf(stdout,"\nend:\t");
	for (i = 0; i < path_n; ++i)
	{
		if (line_end[i]!=0)
			fprintf(stdout, "%d ", line[i][line_end[i]]+hash_len-1);
	}fprintf(stdout, "\noffset:\t");
	for (i = 0; i<path_n; ++i)
	{
		if (line_end[i] != 0)
			fprintf(stdout, "%d ", line[i][0]);
	}
	//fprintf(stdout, "\n");
	for (i = 0; i < path_n; ++i)
		hash_DP_init(h_node, line, line_end, i, hash_len);
	//calculate node value
	for (i = 1; i < path_n; ++i)
		hash_update_line(h_node, line, line_end, i, hash_len);
	//find backtrack node(longest cover-length)
	int max_node_i = -1, max_cover_len = 0;
	for (i = 0; i < path_n; ++i)
	{
		if ((h_node[i].line.ref_cover_len + h_node[i].line.read_cover_len) > max_cover_len)
		{
			max_node_i = i;
			max_cover_len = h_node[i].line.ref_cover_len + h_node[i].line.read_cover_len;
		}
	}
	for (i = 0; i < h_node[max_node_i].line.n_line; ++i)
		copy_flag[h_node[max_node_i].line.line_i[i]] = 1;

	line_end[path_n] = 0;
	for (i = 0; i < path_n; ++i)
	{
		if (copy_flag[i] == 1) hash_copy_line(line, line_end, i, path_n);// XXX len = 3 * line_num
	}

	return path_n;
}*/


/*int hash_copy_main_line(int **line, int *line_end, int m_i, int **mini_line, int *mini_line_end, int mini_i)
{
	if (mini_line_end[mini_i] == 0) return 0;

	int i, j;
	for (i = 0; i < line_end[m_i]; ++i)
	{
		if ((line[m_i][i*3+1] > mini_line[mini_i][1]) || (line[m_i][i*3+1] == mini_line[mini_i][1] && line[m_i][i*3] > mini_line[mini_i][0]))
		{
			break;	
		}
	}
	for (j = line_end[m_i]+mini_line_end[mini_i]-1; j > i+mini_line_end[mini_i]-1; --j)
	{
		line[m_i][j*3] = line[m_i][(j-mini_line_end[mini_i])*3];
		line[m_i][j*3+1] = line[m_i][(j-mini_line_end[mini_i])*3+1];
		line[m_i][j*3+2] = line[m_i][(j-mini_line_end[mini_i])*3+2];
	}
	for (j = i; j < i + mini_line_end[mini_i]; ++j)
	{
		line[m_i][j*3] = mini_line[mini_i][(j-i)*3];
		line[m_i][j*3+1] = mini_line[mini_i][(j-i)*3+1];
		line[m_i][j*3+2] = mini_line[mini_i][(j-i)*3+2];
	}
	line_end[m_i] += mini_line_end[mini_i];
	return 0;
}*/

	//line: i
	//MEM, mis-match allowed seeds combining
	//XXX mis-match have negative effect
	//	int hash_main_line(int *hash_pos, int *start_a, int *len_a, int ref_len, int read_len, int hash_len, int **line, int *line_end, hash_DP_node *h_node)
	//	{
	//		int i, j, k, path_n;
	//		int flag, con_flag;
	//		int path_i[500], tmp, copy_flag[500];
	//		int min_len=len_a[0];
	//
	//		//head hash-node line and tail hash_node line
	//		line[0][0] = 0;
	//		line_end[0] = 0;
	//		line[1][0] = read_len - ref_len;
	//		line_end[1] = 0;
	//		path_n = 2;
	//
	//		//filter out unique-aln seeds
	//		for (i = 0; i <= ref_len-hash_len; ++i) 
	//		{ 
	//			if (len_a[i] == 1)
	//			{
	//				flag = 0;
	//				for (j = path_n-1; j >= 0; --j)
	//				{
	//					hash_main_dis(hash_pos, start_a, i, 0, line[j][0], &con_flag);
	//					if (con_flag == F_MATCH)
	//					{
	//						flag = 1;
	//						line[j][line_end[j]+1] = i;
	//						++line_end[j];
	//						//XXX
	//						/*if (j != (path_n-1))
	//						  {
	//						  for (k = j+1; k != path_n; ++k)
	//						  line_end[k] = 0;
	//						  path_n = j+1;
	//						  }*/
	//						break;
	//					}
	//					//XXX
	//					if (line_end[j] > 5) break;
	//				}
	//				if (flag == 0)
	//				{
	//					hash_new_line(hash_pos[start_a[i]]-i, i, line, line_end, path_n);
	//					++path_n;
	//				}
	//				min_len = 1;
	//			}
	//			else if (min_len != 1 && len_a[i] < min_len && len_a[i] >=2)
	//			{
	//				min_len = len_a[i];	
	//			}
	//		}
	//		//filter out min_len-aln seeds
	//		if (path_n == 2 && line_end[0] == 0 && line_end[1] == 0)
	//		{
	//			//line_end[0] = 0;//no hash-node is uniquely aligned to read
	//			for (i = 0; i <= ref_len-hash_len; ++i)
	//			{
	//				if (len_a[i] == min_len)
	//				{
	//					for (k = 0; k < min_len; ++k)
	//					{
	//						flag = 0;
	//						for (j = path_n-1; j >=0; --j)
	//						{
	//							hash_main_dis(hash_pos, start_a, i, k, line[j][0], &con_flag);
	//							if (con_flag == F_MATCH)
	//							{
	//								flag = 1;
	//								line[j][line_end[j]+1] = i;
	//								++line_end[j];
	//								break;
	//							}
	//							//if (line_end[j] > 5) break;
	//						}
	//						if (flag == 0)
	//						{
	//							hash_new_line(hash_pos[start_a[i]+k]-i, i, line, line_end, path_n);
	//							++path_n;
	//						}
	//					}
	//				}
	//			}
	//			//return 0;
	//		}
	//
	//		//add multi-aln seed to the unique(min_len)-aln seeds lines.
	//		for (i = 0; i <= ref_len-hash_len; ++i)
	//		{
	//			if (len_a[i] > min_len)
	//			{
	//				//XXX k or j ?
	//				for (k = 0; k < path_n; ++k)
	//				{
	//					for (j = 0; j < len_a[i]; ++j)
	//					{
	//						hash_main_dis(hash_pos, start_a, i, j, line[k][0], &con_flag);
	//						if (con_flag == F_MATCH)
	//						{
	//							hash_add_line(line, line_end, k, i);
	//							break;
	//						}
	//					}
	//				}
	//				/*
	//				   for (j = 0; j < len_a[i]; ++j)
	//				   {
	//				   for (k = 0; k < path_n; ++k)
	//				   {//int hash_get_dis(int *hash_pos, int *start_a, int a, int a_i, int b, int b_i, int hash_len, int *con_flag)
	//				//if (i < line[line_end[path_i[k]]-1])
	//				hash_main_dis(hash_pos, start_a, line[path_i[k]][0], i, j, &con_flag);
	//				if (con_flag == F_MATCH)
	//				{
	//				hash_add_line(line, line_end, path_i[k], i);
	//				break;
	//				}
	//				}
	//				}
	//				*/
	//			}
	//		}
	//		hash_cut_line(line, line_end, &path_n, hash_len);
	//
	//		for (i = 0; i < path_n; ++i) { path_i[i] = i; copy_flag[i] = 0; }
	//		//XXX sort by the start postion of every node-line.
	//		//remove the line whose line_end[] = 0
	//		for (i = 0; i < path_n-1; ++i)
	//		{
	//			for (j = i+1; j < path_n; ++j)
	//			{
	//				if (line_end[i] == 0 || (line[j][1] < line[i][1]) || ((line[j][1] == line[i][1]) && (line[j][0] < line[i][0])))
	//					hash_swap_line(line, line_end, i, j);
	//			}
	//		}
	//
	//		/*fprintf(stdout, "all line:\nstart: ");
	//		for (i = 0; i < path_n; ++i)
	//		{
	//			printf("%d, ", line[i][1]);
	//		}
	//		printf("\nend: ");
	//
	//		for (i = 0; i < path_n; ++i)
	//		{
	//			printf("%d, ", line[i][line_end[i]]+hash_len-1);
	//		}
	//		printf("\noffset: ");
	//		for (i = 0; i < path_n; ++i)
	//		{
	//			printf("%d, ", line[i][0]);
	//		}
	//		printf("\n");*/
	//
	//		//DP algorithm
	//		//initialization of DP nodes
	//		//XXX node-lines should be sorted by the start postion
	//		for (i = 0; i < path_n; ++i)
	//			if (line_end[i] != 0) hash_DP_init(h_node, line, line_end, i, hash_len);
	//		//calculate node value
	//		for (i = 1; i < path_n; ++i)
	//			if (line_end[i] != 0) hash_update_line(h_node, line, line_end, i, hash_len);
	//		//find backtrack node(longest cover-length)
	//		int max_node_i = -1, max_cover_len = 0;
	//		for (i = 0; i < path_n; ++i)
	//		{
	//			if (line_end[i] == 0) continue;
	//			if ((h_node[i].line.ref_cover_len + h_node[i].line.read_cover_len) > max_cover_len)
	//			{
	//				max_node_i = i;
	//				max_cover_len = h_node[i].line.ref_cover_len + h_node[i].line.read_cover_len;
	//			}
	//		}
	//		for (i = 0; i < h_node[max_node_i].line.n_line; ++i)
	//			copy_flag[h_node[max_node_i].line.line_i[i]] = 1;
	//		line_end[path_n] = 0;
	//		for (i = 0; i < path_n; ++i)
	//		{
	//			if (copy_flag[i] == 1) hash_copy_line(line, line_end, i, path_n);// XXX len = 3 * line_num
	//		}
	//
	//		//for existing blank, use mini-local_main_line() to fill it
	//		int **_line, *_line_end;
	//		int _ref_start, _ref_end, _read_start, _read_end;
	//		hash_DP_node *_h_node;
	//
	//		_line = (int**)malloc((ref_len-hash_len+1) * sizeof(int*));
	//		_h_node = (hash_DP_node*)malloc((ref_len-hash_len+1) * sizeof(hash_DP_node));
	//		for (i = 0; i <= ref_len-hash_len; ++i) {
	//			_line[i] = (int*)malloc((ref_len-hash_len+1+1) * sizeof(int));
	//			//XXX 100 * sizeof
	//			_h_node[i].line.line_i = (int32_t*)malloc(100 * sizeof(int32_t));
	//		}
	//		_line_end = (int*)malloc((ref_len-hash_len+1) * sizeof(int));
	//
	//		int line_1, line_2, mini_i;
	//
	//		//blank between left bound and head node-line
	//		line_2 = h_node[max_node_i].line.line_i[0];
	//		if (line[line_2][1] > hash_len && (line[line_2][1] + line[line_2][0]) > hash_len)//blank exists XXX > seed_len
	//		{
	//			_ref_start = 0;
	//			_ref_end = line[line_2][1]-1;
	//			_read_start = 0;
	//			_read_end = _ref_end + line[line_2][0];
	//			mini_i = mini_main_line(hash_pos, start_a, len_a, _ref_start, _ref_end, _read_start, _read_end, hash_len, _line, _line_end, _h_node);
	//			if (mini_i != -1 && _line_end[mini_i] != 0)
	//				hash_copy_main_line(line, line_end, path_n, _line, _line_end, mini_i);
	//		}
	//		//blank between node-lines
	//		for (i = 0; i < h_node[max_node_i].line.n_line - 1; ++i)
	//		{
	//			line_1 = h_node[max_node_i].line.line_i[i];
	//			line_2 = h_node[max_node_i].line.line_i[i+1];
	//
	//			if ((line[line_1][line_end[line_1]] + hash_len - 1) < (line[line_2][1]-1-hash_len) //blank exists
	//					&& (line[line_1][line_end[line_1]]+hash_len-1+line[line_1][0]) < (line[line_2][1]+line[line_2][0]-1-hash_len))
	//			{
	//				_ref_start = line[line_1][line_end[line_1]] + hash_len;
	//				_ref_end = line[line_2][1] - 1;
	//				_read_start = _ref_start + line[line_1][0];
	//				_read_end = _ref_end + line[line_2][0];
	//				mini_i = mini_main_line(hash_pos, start_a, len_a, _ref_start, _ref_end, _read_start, _read_end, hash_len, _line, _line_end, _h_node);
	//				if (mini_i != -1 && _line_end[mini_i] != 0)
	//					hash_copy_main_line(line, line_end, path_n, _line, _line_end, mini_i);
	//			}
	//		}
	//		//blank between tail node-line and right bound
	//		line_1 = line_2;
	//		_ref_start = line[line_1][line_end[line_1]] + hash_len;
	//		_ref_end = ref_len - 1;
	//		_read_start = _ref_start + line[line_1][0];
	//		_read_end = read_len - 1;
	//		mini_i = mini_main_line(hash_pos, start_a, len_a, _ref_start, _ref_end, _read_start, _read_end, hash_len, _line, _line_end, _h_node);
	//		if (mini_i != -1 && _line_end[mini_i] != 0)
	//			hash_copy_main_line(line, line_end, path_n, _line, _line_end, mini_i);
	//		//free variable
	//		
	//		for (i = 0; i <= ref_len-hash_len; ++i) { free(_line[i]); free(_h_node[i].line.line_i); }
	//		free(_line); free(_line_end); free(_h_node);
	//
	//		//print filtered node-line
	//		/*fprintf(stdout, "main line: %d\nstart: ", max_node_i);
	//		for (i = 0; i < line_end[path_n]; ++i)
	//		{
	//			printf("%d, ", line[path_n][i*3+1]);
	//		} printf("\nend: ");
	//		for (i = 0; i < line_end[path_n]; ++i)
	//		{
	//			printf("%d, ", line[path_n][i*3+1]+line[path_n][i*3+2]+hash_len - 2);
	//		} printf("\noffset: ");
	//		for (i = 0; i < line_end[path_n]; ++i)
	//		{
	//			printf("%d, ", line[path_n][i*3]);
	//		} printf("\n");*/
	//
	//		return path_n; 
	//	}

/*int hash_add_path(int *hash_pos, int *start_a, int *len_a, 
				  int hash_len, 
				  path_msg **path, int *price_n, 
				  int start, int end, int rev)
{
	if (start == end) return 0;
	else if (start > end) {fprintf(stderr, "[hash_add_path] error: start > end. %d %d\n", start, end); exit(-1);}

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
			for (j = 0; j < len_a[i]; ++j)
			{
				start_dis = hash_get_dis(hash_pos, start_a, anchor, 0, i, j, hash_len, &con_flag);
				if (con_flag != F_UNCONNECT) {
					last_dis = hash_get_dis(hash_pos, start_a, last_i, last_j, i, j, hash_len, &last_flag);
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
							hash_get_dis(hash_pos, start_a, k, l, i, tmp_j, hash_len, &last_flag);
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
			for (j = 0; j < len_a[i]; ++j)
			{
				start_dis = hash_get_dis(hash_pos, start_a, start, 0, i, j, hash_len, &con_flag);
				if (con_flag != F_UNCONNECT) 
				{
					hash_get_dis(hash_pos, start_a, i, j, end, 0, hash_len, &end_flag);
					if (end_flag != F_UNCONNECT) {
						last_dis = hash_get_dis(hash_pos, start_a, last_i, last_j, i, j, hash_len, &last_flag);
						//XXX aln-res are both F_MATCH, then compare on the cigar_len and edit-dis.
						//  if (last_flag == F_MATCH) {
						//  path[i][j].from.x = last_i;
						//  path[i][j].from.y = last_j;
						//  path[i][j].flag = F_MATCH;
						//  last_i = i; last_j = j; last_start_dis = start_dis;
						//  price_n[i] = j+1;
						//  flag = 0;
						//  break;
						//  }
						//  else 
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
							hash_get_dis(hash_pos, start_a, k, l, i, tmp_j, hash_len, &last_flag);
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
		else {fprintf(stderr, "[hash_add_path] Bug: rev error.\n"); exit(-1);}
	}
	return 1;
}*/

/*int hash_frag_set(hash_frag *h_frag, int ref_bound, int read_bound, int flag, int hash_len)
{
	if (flag == HASH_FRAG_START)
	{
		h_frag->hfa[h_frag->frag_num].ref_start = ref_bound;
		h_frag->hfa[h_frag->frag_num].read_start = read_bound;
		++h_frag->frag_num;
	}
	else if (flag == HASH_FRAG_END)
	{
		h_frag->hfa[h_frag->frag_num].ref_end = ref_bound + hash_len -1;
		h_frag->hfa[h_frag->frag_num].read_end = read_bound + hash_len - 1;
	}
	else	//HASH_FRAG_SEED
	  {

	  }
	return 0;
}*/

/*int hash_backtrack(path_msg **path, int *price_n, int *hash_pos, int *start_a, int n_hnode, int hash_len, hash_frag *h_frag)
{
	//Determin the start point of backtrack.
	int i;//(path[i][min_i])
	for (i = n_hnode-1; i >= 0; i--)
		if (price_n[i] > 0)	break;

	//backtrack from (i, min_i)
	int last_x = i, last_y = price_n[i]-1, frag_num = 0, tmp, flag = 0;
	int end_x, start_x;
	//first end
	//frag_set_msg(a_msg, last_x, last_y, 1, f_msg, frag_num, seed_len);
	hash_frag_set(h_frag, last_x, hash_pos[start_a[last_x]+last_y], HASH_FRAG_END, hash_len);
	end_x = last_x;
	//fprintf(stdout, "    end     ref %d read %d\n", last_x+9, hash_pos[start_a[last_x]+last_y]+9);
	while (path[last_x][last_y].flag != F_PATH_END) {
		if (path[last_x][last_y].flag != F_MATCH) {
			//start
			start_x = last_x;
			if (start_x == end_x) {
				//fprintf(stdout, "    flase start ref %d read %d\n", last_x, hash_pos[start_a[last_x]+last_y]);
			}
			else {
				//frag_set_msg(a_msg, last_x, last_y, 0, f_msg, frag_num, seed_len);
				hash_frag_set(h_frag, last_x, hash_pos[start_a[last_x]+last_y], HASH_FRAG_START, hash_len);
				//fprintf(stdout, "    start   ref %d read %d\n", last_x, hash_pos[start_a[last_x]+last_y]);
				++frag_num;
			}
			flag = 1;
		}
		tmp = last_x;
		last_x = path[last_x][last_y].from.x;
		last_y = path[tmp][last_y].from.y;

		if (flag == 1) {
			//next end
			//frag_set_msg(a_msg, last_x, last_y, 1, f_msg, frag_num, seed_len);
			hash_frag_set(h_frag, last_x, hash_pos[start_a[last_x]+last_y], HASH_FRAG_END, hash_len);
			end_x = last_x;

			//fprintf(stdout, "    end     ref %d read %d\n", last_x+9, hash_pos[start_a[last_x]+last_y]+9);
			flag = 0; 
		} else {;
			//frag_set_msg(a_msg, last_x, last_y, 2, f_msg, frag_num, seed_len);	//seed
			hash_frag_set(h_frag, last_x, hash_pos[start_a[last_x]+last_y], HASH_FRAG_SEED, hash_len);
			//printf("seed: ref %d read %d\n", last_x, hash_pos[start_a[last_x]+last_y]);
		}
	}
	//start
	//frag_set_msg(a_msg, last_x, last_y, 0, f_msg, frag_num, seed_len);
	hash_frag_set(h_frag, last_x, hash_pos[start_a[last_x]+last_y], HASH_FRAG_START, hash_len);
	//fprintf(stdout, "    start   ref %d read %d\n", last_x, hash_pos[start_a[last_x]+last_y]);
	return 1;
}*/

//return the overlap len of two nodes 
//return value >= 0
int make_indel_cigar(int ref_left, int read_left, int ref_right, int read_right, int *clen, uint32_t **cigar)
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
	}
	else if (len < 0)
	{
		(*clen) = 1;
		(*cigar)[0] = ((0-len) << 4) + CINS;
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
		if (_head) 
		{
			h_node[0][0] = (hash_dp_node){{-1,0}, 0-hash_len, 0, 0, 0, F_MATCH, MIN_FLAG};
			if (_tail) { h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, ref_len, read_len-ref_len, 0, 0, F_UNMATCH, MIN_FLAG}; bound_flag=1;} //split_offset=read_len-ref_len; }
			else { h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG}; }
		}
		else
		{
			h_node[0][0] = (hash_dp_node){{-1,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG};
			if (_tail) { h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, ref_len, read_len-ref_len, 0, 0, F_MATCH, MIN_FLAG}; }// split_offset = read_len-ref_len; } 
			else h_node[hash_seed_n+1][0] = (hash_dp_node){{0,0}, -1, -1, 0, 0, F_MATCH, UNLIMITED_FLAG};
		}

		//seed init
		for (i = 1; i <= hash_seed_n; ++i)
		{
			if (len_a[i] == min_len)
			{
				hash_dp_init(h_node, hash_pos, start_a, len_a, i, (i-1)*hash_step/*ref offset*/, head, hash_len, hash_step, MIN_FLAG, split_offset, 1/*bound_flag*/, bound_len);
				min_exist = 1;
			}
			else hash_dp_init(h_node, hash_pos, start_a, len_a, i, (i-1)*hash_step/*ref offset*/, head, hash_len, hash_step, MULTI_FLAG, split_offset, 1/*bound_flag*/, bound_len);
		}
	}
	//dp update and backtrack
	{
		//min update and backtrack
		if (min_exist)
		{
			//min extend
#ifndef __NEW__
			for (i = 1; i <= hash_seed_n; ++i)
			{
				if (len_a[i] > min_len)
					hash_min_extend(h_node, len_a, i, hash_seed_n+2, min_len, MIN_FLAG, hash_len, hash_step, split_offset);
			}
#endif
#ifdef __NEW__
			if (_head) hash_min_extend(h_node, len_a, head.x, head.y, hash_seed_n+2, min_len, MIN_FLAG, hash_len, hash_step, split_offset);
			for (i = 1; i <= hash_seed_n; ++i)
			{
				if (len_a[i] <= min_len && len_a[i] > 0)
				{
					for (j = 0; j < len_a[i]; ++j)
					{
						hash_min_extend(h_node, len_a, i, j, hash_seed_n+2, min_len, MIN_FLAG, hash_len, hash_step, split_offset);
					}
				}
			}
			if (_tail) hash_min_extend(h_node, len_a, tail.x, tail.y, hash_seed_n+2, min_len, MIN_FLAG, split_offset, hash_len, hash_step);
#endif
			//min update
			for (i = 2; i <= hash_seed_n; ++i)
			{
				for (j = 0; j < len_a[i]; ++j)
				{
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
				if (left.x == head.x)
					break;
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
			//DEBUG
                /*if (node_i == 0)
                {
                    printf("min-update\n");
                    for (i = 0; i <= hash_seed_n+1; ++i)
                    {
                        for (j = 0; j < len_a[i]; ++j)
                        {
                            printf("%d %d: from %d %d ref_i %d offset %d node_n %d score %d dp_flag %d\n",i,j, h_node[i][j].from.x, h_node[i][j].from.y, h_node[i][j].ref_i, h_node[i][j].offset, h_node[i][j].node_n, h_node[i][j].score, h_node[i][j].dp_flag);
                        }
                    }
                }*/
				/*{
					printf("min-update\n");
					for (i = 0; i < node_i; ++i)
						printf("%d: from(%d %d) ref_i %d offset %d score %d node_n %d match_f %d dp_f %d\n", i, h_node[line[i].x][line[i].y].from.x,h_node[line[i].x][line[i].y].from.y, h_node[line[i].x][line[i].y].ref_i, h_node[line[i].x][line[i].y].offset, h_node[line[i].x][line[i].y].score, h_node[line[i].x][line[i].y].node_n, h_node[line[i].x][line[i].y].match_flag, h_node[line[i].x][line[i].y].dp_flag);
				}*/

			return node_i;
		}
		//whole-multi update
		else	
		{
			for (i = 2; i <= hash_seed_n; ++i)
			{
				for (j = 0; j < len_a[i]; ++j)
				{
					if (h_node[i][j].dp_flag == MULTI_FLAG)
						hash_dp_update(h_node, len_a, i, j, 1/*update start pos*/, hash_len, hash_step, MULTI_FLAG, split_offset, bound_flag, bound_len);
				}
			}
			hash_dp_update(h_node, len_a, tail.x, tail.y, 1/*update start pos*/, hash_len, hash_step, MULTI_FLAG, split_offset, bound_flag, bound_len);
			node_i = h_node[tail.x][tail.y].node_n-1;
				//DEBUG
				/*if (node_i == -1)
				{
                    printf("whole-update\n");
                    for (i = 0; i <= hash_seed_n+1; ++i)
                    {
                        for (j = 0; j < len_a[i]; ++j)
                        {
                            printf("%d %d: from %d %d ref_i %d offset %d node_n %d score %d dp_flag %d\n",i,j, h_node[i][j].from.x, h_node[i][j].from.y, h_node[i][j].ref_i, h_node[i][j].offset, h_node[i][j].node_n, h_node[i][j].score, h_node[i][j].dp_flag);
                        }
                    }
                }*/

			int last_x, last_y;
			//right = tail;
			i = h_node[tail.x][tail.y].from.x;
			j = h_node[tail.x][tail.y].from.y;
			//left = h_node[right.x][right.y].from;
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
				//DEBUG
				/*{
					printf("whole-update\n");
					for (i = 0; i < h_node[tail.x][tail.y].node_n; ++i)
						printf("%d: from(%d %d) ref_i %d offset %d score %d node_n %d match_f %d dp_f %d\n", i, h_node[line[i].x][line[i].y].from.x,h_node[line[i].x][line[i].y].from.y, h_node[line[i].x][line[i].y].ref_i, h_node[line[i].x][line[i].y].offset, h_node[line[i].x][line[i].y].score, h_node[line[i].x][line[i].y].node_n, h_node[line[i].x][line[i].y].match_flag, h_node[line[i].x][line[i].y].dp_flag);
				}*/
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
int hash_split_map(uint32_t **split_cigar, uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, int hash_len, int key_len, uint32_t *hash_num, uint64_t **hash_node, int32_t *hash_pos)
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
int hash_split_map(uint32_t **split_cigar, int *split_clen, int *split_m,
				   uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
				   int hash_len, int hash_step, int key_len, 
				   uint32_t *hash_num, uint64_t **hash_node, int32_t *hash_pos, 
				   int _head, int _tail)
{
	int i, res=0;
	uint8_t *ref_query;
	int q_key_int, q_kmer_int, node_i;

    //un-overlap XXX
    int hash_seed_n = (ref_len-hash_len)/hash_step + 1;
	int *start_a = (int*)malloc((hash_seed_n + 2) * sizeof(int));	//2: for head and tail
	int *len_a = (int*)malloc((hash_seed_n + 2) * sizeof(int));
	if (start_a == NULL || len_a == NULL) {fprintf(stderr, "[hash_split_map] Not enougy memory.(ref_len %d, read_len %d)\n", ref_len, read_len); exit(-1);}

	(*split_clen) = 0;

	//hash map, restore hash map result
	{
		len_a[0] = 1;		//for head node
		for (i = 0; i < ref_len; ++i) if (ref_seq[i] > 3) {free(start_a); free(len_a); return 0;}
		for (i = 0; i <= ref_len - hash_len; i+=hash_step)
		{
			ref_query = ref_seq + i;
			if (hash_calcu(&q_key_int, &q_kmer_int, ref_query, hash_len, key_len) == 1)	// hash 'N'
			{
				start_a[i/hash_step+1] = 0;
				len_a[i/hash_step+1] = 0;
			}
			else
			{
				if (hash_hit(hash_num, hash_node, &node_i, q_key_int, q_kmer_int, hash_len-key_len) == 1)
				{
					start_a[i/hash_step + 1] = ((hash_node[q_key_int][node_i] & 0xffff0000) >> 16);	//for hit seeds
					len_a[i/hash_step + 1] = (hash_node[q_key_int][node_i] & 0xffff);
				}
				else len_a[i/hash_step + 1] = 0;	//for un-hit seeds
			}
		}
		len_a[i/hash_step + 1] = 1;	//for tail node
	}
	//hash-node DP, select and arrange seeds's order on ref
		line_node *line;
		hash_dp_node **h_node;
		int m_len;
	{
		line = (line_node*)malloc(hash_seed_n * sizeof(line_node));	//store the path of DP result
		h_node = (hash_dp_node **)malloc((hash_seed_n + 2) * sizeof(hash_dp_node*));
		for (i = 0; i < hash_seed_n+2; ++i)
			h_node[i] = (hash_dp_node*)malloc(len_a[i] * sizeof(hash_dp_node));
		if (line == NULL || h_node == NULL) {fprintf(stderr, "[hash_split_map] Not enougy memory.\n"); exit(-1);} 

		m_len = hash_main_line(hash_pos, start_a, len_a, ref_len, read_len, hash_seed_n, hash_len, hash_step, h_node, line, _head, _tail);
		/*if (m_len == 0) 
		{
			for (i = 0; i < ref_len; ++i) printf("%c", "ACTG"[ref_seq[i]]);
			fprintf(stdout, "\n");
			for (i = 0; i < read_len; ++i) printf("%c", "ACTG"[read_seq[i]]);
			fprintf(stdout, "\n[hash_split_map] error.\n"); 
			goto _FREE_V;
		}*/
	}

	int _q_len, _t_len, _clen, _b_w, _score; 
	uint32_t *_cigar=0;
	if (m_len > 0)
	{
		//fill blank with generated SV and SW, return the whole cigar
			uint32_t *g_cigar;
			g_cigar = (uint32_t*)malloc(100 * sizeof(uint32_t));
			int tail_in, head_in;

			//XXX hash read len == 0 XXX
		//1. fix the region between left bound and first line
			if (_head)
			{
				if (h_node[line[0].x][line[0].y].ref_i != 0 && (h_node[line[0].x][line[0].y].ref_i + h_node[line[0].x][line[0].y].offset) != 0)	//blank exists
				{
					/*_t_len = h_node[line[0].x][line[0].y].ref_i;
					_q_len = h_node[line[0].x][line[0].y].ref_i + h_node[line[0].x][line[0].y].offset;
					_b_w = abs(_t_len - _q_len);
					_score = ksw_global(_q_len, read_seq, _t_len, ref_seq, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);*/

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
					if (_clen/2 > 10 && _score < 0 && _q_len > 100 && _t_len > 100)
					{
						res = 1;
						_cigar[0] = CSOFT_CLIP;
						_cigar[1] = (_q_len << 4) | CINS;
						_cigar[2] = (_t_len << 4) | CDEL;
						_cigar[3] = CSOFT_CLIP;
						_clen = 4;
					}
					_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
					/*if ((*split_cigar)[*split_clen-1] & 0xf != CMATCH) tail_out = 0;
					else
					{
						tail_out = ((*split_cigar)[*split_clen-1] >> 4);
						if (tail_out > 10) tail_out = 10;
					}*/
					//_cigar[0] = ((-tail_out) << 4 | CMATCH); _clen = 1;
					//_push_cigar(split_cigar, split_clen, _cigar, _clen);
					free(_cigar);
				}
				else//no blank, add SV cigar
				{	
					//         ref_left, read_left, ref_right,           read_right
					make_indel_cigar(-1, -1,  h_node[line[0].x][line[0].y].ref_i, h_node[line[0].x][line[0].y].ref_i+h_node[line[0].x][line[0].y].offset, &_clen, &g_cigar);
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
						if (g_cigar[0] >>4 > 10000) { fprintf(stderr, "match both tail:%d overlap: %d\t", tail_in, overlap); printcigar(g_cigar, _clen); printf("\n"); }

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
						/*if (_score < 0)
						{
							int j;
							printf("hash old ref:\t"); for (j = 0; j < _t_len; ++j) printf("%c", "ACGT"[(ref_seq+h_node[line[i].x][line[i].y].ref_i+hash_len-head_in)[j]]); printf("\n");
							printf("hash old read:\t"); for(j = 0; j < _q_len; ++j) printf("%c", "ACGT"[(read_seq+h_node[line[i].x][line[i].y].ref_i+hash_len+h_node[line[i].x][line[i].y].offset-head_in)[j]]); printf("\n");
							printf("hash old both cigar:\t");
							printcigar(_cigar, _clen);
							k = 0;
							for (j = 0; j < _clen; ++j)
							{
								if ((_cigar[j] & 0xf) != CMATCH)
									k++;
							}
							fprintf(stdout, "\t%d\t%d\n", k, _clen);

							//	_t_len = h_node[line[i+1].x][line[i+1].y].ref_i - h_node[line[i].x][line[i].y].ref_i - hash_len + 2 * hash_len;
							//	_q_len = _t_len + h_node[line[i+1].x][line[i+1].y].offset - h_node[line[i].x][line[i].y].offset;
							//	_b_w = abs(_t_len - _q_len);
							//	_score = ksw_global(_q_len, read_seq+h_node[line[i].x][line[i].y].ref_i+hash_len+h_node[line[i].x][line[i].y].offset-hash_len, _t_len, ref_seq+h_node[line[i].x][line[i].y].ref_i+hash_len-hash_len, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
							//	printf("both hash new ref:\t"); for (j = 0; j < _t_len; ++j) printf("%c", "ACGT"[(ref_seq+h_node[line[i].x][line[i].y].ref_i+hash_len-hash_len)[j]]); printf("\n");
							//	printf("hash new read:\t"); for(j = 0; j < _q_len; ++j) printf("%c", "ACGT"[(read_seq+h_node[line[i].x][line[i].y].ref_i+hash_len+h_node[line[i].x][line[i].y].offset-hash_len)[j]]); printf("\n");
							//	printf("hash new cigar:\t");
							//	printcigar(_cigar, _clen);
							//	fprintf(stdout, "\n");
						}*/
						if (_clen/2 > 10 && _score < 0 && _q_len > 100 && _t_len > 100)
						{
							res = 1;
							_cigar[0] = CSOFT_CLIP;
							_cigar[1] = (_q_len << 4) | CINS;
							_cigar[2] = (_t_len << 4) | CDEL;
							_cigar[3] = CSOFT_CLIP;
							_clen = 4;
						}
						_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
						overlap = 0;
						free(_cigar);
					}
					else//no blank, add SV cigar
					{
						tail_in = 0;
						overlap = make_indel_cigar(h_node[line[i].x][line[i].y].ref_i+hash_len-1, h_node[line[i].x][line[i].y].ref_i+hash_len-1+h_node[line[i].x][line[i].y].offset, h_node[line[i+1].x][line[i+1].y].ref_i, h_node[line[i+1].x][line[i+1].y].ref_i+h_node[line[i+1].x][line[i+1].y].offset, &_clen, &g_cigar);
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
					/*if (_score < 0)
					{
						int j;
						printf("hash old ref:\t"); for (j = 0; j < _t_len; ++j) printf("%c", "ACGT"[(ref_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len - head_in)[j]]); printf("\n");
						printf("hash old read:\t");for (j = 0; j < _q_len; ++j) printf("%c", "ACGT"[(read_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+h_node[line[m_len-1].x][line[m_len-1].y].offset+hash_len - head_in)[j]]); printf("\n");
						printf("hash old tail cigar:\t");
						printcigar(_cigar, _clen);
						k = 0;
						for (j = 0; j < _clen; ++j)
						{
							if ((_cigar[j] & 0xf) != CMATCH)
								k++;
						}
						fprintf(stdout, "\t%d\t%d\n", k, _clen);
						//	_t_len = ref_len - h_node[line[m_len-1].x][line[m_len-1].y].ref_i - hash_len+hash_len;
						//	_q_len = read_len - h_node[line[m_len-1].x][line[m_len-1].y].ref_i - hash_len - h_node[line[m_len-1].x][line[m_len-1].y].offset + hash_len;
						//	_b_w = abs(_t_len - _q_len);
						//	_score = ksw_global(_q_len, read_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+h_node[line[m_len-1].x][line[m_len-1].y].offset+hash_len-hash_len, _t_len, ref_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len-hash_len, 5, bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
						//	printf("tail hash new ref:\t"); for (j = 0; j < _t_len; ++j) printf("%c", "ACGT"[(ref_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len-hash_len)[j]]); printf("\n");
						//	printf("hash new read:\t");for (j = 0; j < _q_len; ++j) printf("%c", "ACGT"[(read_seq+h_node[line[m_len-1].x][line[m_len-1].y].ref_i+h_node[line[m_len-1].x][line[m_len-1].y].offset+hash_len-hash_len)[j]]); printf("\n");
						//	printf("hash new cigar:\t");
						//	printcigar(_cigar, _clen);
						//	fprintf(stdout, "\n");
					}*/

					if (_clen/2 > 10 && _score < 0 && _q_len > 100 && _t_len > 100)
					{
						res = 1;
						_cigar[0] = CSOFT_CLIP;
						_cigar[1] = (_q_len << 4) | CINS;
						_cigar[2] = (_t_len << 4) | CDEL;
						_cigar[3] = CSOFT_CLIP;
						_clen = 4;
					}
					_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
					free(_cigar);
				}
				else //no blank, add SV cigar
				{
					make_indel_cigar(h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len-1, h_node[line[m_len-1].x][line[m_len-1].y].ref_i+hash_len-1+h_node[line[m_len-1].x][line[m_len-1].y].offset, ref_len, read_len, &_clen, &g_cigar);
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
			if (_clen/2 > 10 && _score < 0 && _q_len > 100 && _t_len > 100)
			{
				res = 1;
				_cigar[0] = CSOFT_CLIP;
				_cigar[1] = (_q_len << 4) | CINS;
				_cigar[2] = (_t_len << 4) | CDEL;
				_cigar[3] = CSOFT_CLIP;
				_clen = 4;
			}
			_push_cigar(split_cigar, split_clen, split_m, _cigar, _clen);
			free(_cigar);
		}
	}
//_FREE_V:
	for (i = 0; i < hash_seed_n+2; ++i) free(h_node[i]);
	free(h_node); free(line); 
	free(start_a); free(len_a);
	return res;
}

//XXX for INS: the ref and read bound maybe the right bound of INS seqs, cause of the INS is just a repeat of seqs nearby.
/*int find_bound(uint8_t *ref_seq, int ref_len, int read_len, int hash_len, int key_len, uint32_t *hash_num, uint64_t **hash_node, int32_t *hash_pos, int *ref_b1, int *ref_b2, int *read_b1, int *read_b2)
{
	int i, res, read_i, un_hit;
	uint8_t *ref_query;

	*read_b1 = *ref_b1 = -1;
	*read_b2 = read_len;
	*ref_b2 = ref_len;
	for (i=0, un_hit=0; i<=ref_len-hash_len && un_hit<hash_len+1; ++i)
	{
		ref_query = ref_seq+i;
		//printf("i: %d ", i);
		res = hash_map(hash_num, hash_node, hash_pos, &read_i, ref_query, i, hash_len, key_len);
		if (res == 2) ++un_hit;
		else if (res == 1 && read_i > (*read_b1+1-hash_len))
		{
			un_hit=0;
			*ref_b1 = i + hash_len - 1;
			*read_b1 = read_i + hash_len - 1;
		}
	}
	for (i=ref_len-hash_len, un_hit=0; i>*ref_b1 && un_hit<hash_len+1; --i)
	{
		ref_query = ref_seq+i;
		//printf("I: %d ", i+read_len-ref_len);
		res = hash_map(hash_num, hash_node, hash_pos, &read_i, ref_query, i+read_len-ref_len, hash_len, key_len);
		if (res == 1 && read_i < *read_b2)
		{
			un_hit = 0;
			if (read_i <= *read_b1) break;
			*ref_b2 = i;
			*read_b2 = read_i;
		}
		else ++un_hit;
	}
	return 0;
}*/

/*int separate_cigar(uint32_t *cigar, int cigar_len, int slen, uint32_t *c1, int *c1_len, uint32_t *c2, int *c2_len)
{
	int i, j, llen = slen;

	*c1_len = *c2_len = 0;
	for (i = 0; i < cigar_len; ++i)
	{
		if ((cigar[i] >> 4) == llen)
		{
			c1[*c1_len] = CIGAR_GEN(llen, cigar[i]&0xf);//(llen << 4) + (int)(cigar[i] & 0xf);
			(*c1_len)++;
			// EQUAL
			*c2_len = cigar_len - *c1_len;
			for (j = 0; j < *c2_len; ++j)
			{
				c2[j] = cigar[i+j+1];
			}
			return 0;
		}
		else if ((cigar[i] >> 4) > llen)
		{
			c1[*c1_len] = CIGAR_GEN(llen, cigar[i]&0xf);//(llen << 4) | (int)(cigar[i] & 0xf);
			(*c1_len)++;
			//NOT EQUAL
			*c2_len = cigar_len - *c1_len + 1;
			c2[0] = CIGAR_GEN((cigar[i]>>4)-llen, cigar[i]&0xf);//(int)(((cigar[i]>>4) - llen) << 4) + (int)(cigar[i] & 0xf);
			for (j = 1; j < *c2_len; ++j)
				c2[j] = cigar[i+j];
			return 0;
		}
		else
		{
			c1[*c1_len] = cigar[i];
			(*c1_len)++;
			llen = llen - (cigar[i]>>4);
		}
	}
	return 0;
}*/

/*int amend_bound_core(uint32_t *c1, int c1_len, uint32_t *c2, int c2_len, 
					 int *sseq_b1, int *sseq_b2, int *lseq_b1, int *lseq_b2)
{
	int i,j;
	uint32_t c;
	for (i = c1_len-1; i>=0; --i)
	{
		if ((c1[i] & 0xf) == CMATCH)
		{
			*sseq_b1 += c1[i]>>4;
			*lseq_b1 += c1[i]>>4;
			for (j = i-1; j >= 0; --j)
			{
				c = c1[j] & 0xf;
				if (c == CMATCH)
				{
					*sseq_b1 += c1[j]>>4;
					*lseq_b1 += c1[j]>>4;
				}
				else if (c == CINS)
					*lseq_b1 += c1[j]>>4;
				else *sseq_b1 += c1[j]>>4;
			}
		}
	}
	for (i = 0; i < c2_len; ++i)
	{
		if ((c2[i] & 0xf) == CMATCH)
		{
			*sseq_b2 -= c2[i]>>4;
			*lseq_b2 -= c2[i]>>4;
			for (j = i+1; j < c2_len; ++j)
			{
				c = c2[j] & 0xf;
				if (c == CMATCH)
				{
					*sseq_b2 -= c2[j]>>4;
					*lseq_b2 -= c2[j]>>4;
				}
				else if (c == CINS)
					*sseq_b2 -= c2[j]>>4;
				else *lseq_b2 -= c2[j]>>4;
			}
		}
	}
	return 0;
}*/

/*int amend_bound(uint8_t *sseq, int *sseq_b1, int *sseq_b2, uint8_t *lseq, int *lseq_b1, int *lseq_b2)
{
	int tlen = *sseq_b2-*sseq_b1-1, qlen,i, b_w, score;
	uint8_t *target = sseq+*sseq_b1+1, *query;
	uint32_t *cigar=0, *c1, *c2;
	int cigar_len, c1_len, c2_len;
	if ((*lseq_b2-*lseq_b1-1) >= 2 * tlen)
	{
		qlen = 2 * tlen;
		query = (uint8_t*)malloc(qlen * sizeof(uint8_t));
		for (i = 0; i < qlen; ++i)
		{
			if (i < tlen) query[i] = lseq[*lseq_b1+1+i];
			else query[i] = lseq[*lseq_b2-(qlen-i)];
		}
		b_w = tlen;
	}
	else
	{
		qlen = *lseq_b2-*lseq_b1-1;
		query = lseq+*lseq_b1+1;
		b_w = qlen - tlen;
	}
	score = ksw_global(qlen, query, tlen, target, 5, bwasw_sc_mat, 5, 2, b_w, &cigar_len, &cigar); 
	//printcigar(cigar, cigar_len);
	c1 = (uint32_t*)malloc(cigar_len * sizeof(uint32_t)); c2 = (uint32_t*)malloc(cigar_len * sizeof(uint32_t));
	separate_cigar(cigar, cigar_len, tlen, c1, &c1_len, c2, &c2_len);

	//printcigar(c1, c1_len);	printcigar(c2, c2_len);
	amend_bound_core(c1, c1_len, c2, c2_len, sseq_b1, sseq_b2, lseq_b1, lseq_b2);

	if (qlen == 2 * tlen) free(query);
	free(cigar); free(c1); free(c2);
	return 0;
}*/

/*int get_res_cigar(uint8_t *sseq, int slen, int sseq_b1, int sseq_b2, 
				  uint8_t *lseq, int llen, int lseq_b1, int lseq_b2, 
				  uint32_t **res_cigar, int *cigar_len, int FLAG)
{
	int8_t cin;
	if (FLAG == 1) cin = CINS; else cin = CDEL;

	int tlen, qlen, b_w, score, c_len, c1_len, c2_len, i;
	uint32_t *c1, *c2;
	tlen = sseq_b1+1; qlen = lseq_b1+1; b_w = abs(tlen-qlen);
	score = ksw_global(qlen, lseq, tlen, sseq, 5, bwasw_sc_mat, 5, 2, b_w, &c1_len, &c1); //printcigar(c1, c1_len);
	tlen = slen - sseq_b2;	qlen = llen - lseq_b2; b_w = abs(tlen -qlen);

	score = ksw_global(qlen, lseq+lseq_b2, tlen, sseq+sseq_b2, 5, bwasw_sc_mat, 5, 2, b_w, &c2_len, &c2); //printcigar(c2, c2_len);

	// *res_cigar = (uint32_t*)malloc((c1_len+c2_len+2)*sizeof(uint32_t));
	if (FLAG == 1) for (c_len = 0; c_len < c1_len; ++c_len) (*res_cigar)[c_len] = c1[c_len];
	else
	{
		for (c_len = 0; c_len < c1_len; ++c_len)
		{
			if ((c1[c_len]&0xf) == CMATCH) (*res_cigar)[c_len] = c1[c_len];
			else if ((c1[c_len]&0xf) == CINS) (*res_cigar)[c_len] = c1[c_len]+1;
			else (*res_cigar)[c_len] = c1[c_len]-1;
		}
	}
	if (((*res_cigar)[c_len-1]&0xf) == CMATCH)
	{
		if (sseq_b2-sseq_b1-1 > 0) (*res_cigar)[c_len-1] |= ((sseq_b2-sseq_b1-1)<<4);
		if (lseq_b2-lseq_b1-sseq_b2+sseq_b1 > 0) (*res_cigar)[c_len++] = CIGAR_GEN(lseq_b2-lseq_b1-sseq_b2+sseq_b1, cin);
	}
	else
	{
		if (sseq_b2-sseq_b1-1 > 0) (*res_cigar)[c_len++] = CIGAR_GEN(sseq_b2-sseq_b1-1, CMATCH);
		if (lseq_b2-lseq_b1-sseq_b2+sseq_b1 > 0)(*res_cigar)[c_len++] = CIGAR_GEN(lseq_b2-lseq_b1-sseq_b2+sseq_b1, cin);
	}
	if (FLAG == 1) for (i = 0; i < c2_len; ++i, ++c_len) (*res_cigar)[c_len] = c2[i];
	else
	{
		for (i = 0; i < c2_len; ++i, ++c_len)
		{
			if ((c2[i]&0xf) == CMATCH) (*res_cigar)[c_len] = c2[i];
			else if ((c2[i]&0xf) == CINS) (*res_cigar)[c_len] = c2[i]+1;
			else (*res_cigar)[c_len] = c2[i]-1;
		}
	}
	*cigar_len = c_len;

	free(c1); free(c2);
	return 0;
}*/

//srand+/- XXX
//main function of split-mapping
//read_seq : seq of read part that contain a breakpoint or more. / read_len
//ref1/2_offset : offset of two end's reference(left most)
//hash_num/hash_node : hash table
//ha_msg : struct for hash-map result
//key_len : length of first path of hash-sequence;/hash_size : the whole number of posibility of key_len bps;
//return: XXX ?

//XXX hash_node: use struct?
//insertion length = read_len - ref_len
int split_insert_map(uint32_t **res_cigar, int *res_len, int *res_m,
					 uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, 
					 int64_t ref_offset, 
					 int hash_len, int hash_step, 
					 uint32_t **hash_num, uint64_t ***hash_node, 
					 int key_len, int hash_size)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	/*
		for (i = 0; i < read_len; ++i) fprintf(stdout, "%d", read_seq[i]);
		fprintf(stdout, "\n");
		for (i = 0; i < ref_len; ++i) fprintf(stdout, "%d", ref_seq[i]);
		fprintf(stdout,"\n");*/
	//init of hash table
	init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_pos, key_len, hash_size);

	int res;
	//for multi-breakpoints
	res = hash_split_map(res_cigar, res_len, res_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, key_len, *hash_num, *hash_node, hash_pos, 1, 1);
	//printcigar(*res_cigar, *res_len); fprintf(stdout,"\n");
	
    //printf("hash ins:\t"); printcigar(*res_cigar, *res_len); printf("\n");
	free(hash_pos);
	return res;
}

//deletion length = ref_len - read_len
int split_delete_map(uint32_t **res_cigar, int *res_len, int *res_m,
					 uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, 
					 int64_t ref_offset, 
					 int hash_len, int hash_step, 
					 uint32_t **hash_num, uint64_t ***hash_node, 
					 int key_len, int hash_size)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	/* 
		for (i = 0; i < read_len; ++i) fprintf(stdout, %d", read_seq[i]);
		fprintf(stdout, "\n");
		for (i = 0; i < ref_len; ++i) fprintf(stdout, "%d", ref_seq[i]);
		fprintf(stdout, "\n");*/
	//init of hash table
	init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_pos, key_len, hash_size);

	int res;
	//for multi-breakpoints
	res = hash_split_map(res_cigar, res_len, res_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, key_len, *hash_num, *hash_node, hash_pos, 1, 1);
	/*printcigar(*res_cigar, *res_len); fprintf(stdout,"\n");*/

	//printf("hash del:\t"); printcigar(*res_cigar, *res_len); printf("\n");
	free(hash_pos);
	return res;
}

int hash_right_bound_map(uint32_t **cigar, int *cigar_len, int *cigar_m,
						 uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
						 uint32_t **hash_num, uint64_t ***hash_node, 
						 int hash_len, int hash_key, int hash_step)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	int hash_size = pow(NT_N, hash_key);
	if (init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_pos, hash_key, hash_size) != 0)
	{
		(*cigar_len) = 0;
		free(hash_pos);
		return 1;
	}
	int res;
	res = hash_split_map(cigar, cigar_len, cigar_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, hash_key, *hash_num, *hash_node, hash_pos, 1, 0);
	int i, readINcigar=0, refINcigar=0;
	for (i = 0; i < (*cigar_len); ++i)
	{
		if (((*cigar)[i] & 0xf) == CMATCH || ((*cigar)[i] & 0xf) == CSOFT_CLIP) readINcigar += ((*cigar)[i] >> 4), refINcigar += ((*cigar)[i] >> 4);
		else if (((*cigar)[i] & 0xf) == CINS) readINcigar += ((*cigar)[i] >> 4);
		else if (((*cigar)[i] & 0xf) == CDEL) refINcigar += ((*cigar)[i] >> 4);
		else 
        {
            fprintf(stderr, "[hash_right_bound_map] cigar error: ");
            printcigar(*cigar, *cigar_len);
            fprintf(stderr, "\n");
            exit(-1);
        }
	}
	if (readINcigar < read_len)     
	{
		int qlen = read_len-readINcigar;
		int tlen = read_len -readINcigar + ((hash_len+read_len-readINcigar)*bwasw_sc_mat[0]-5-1)/2;
		tlen = tlen < (ref_len - refINcigar) ? tlen : (ref_len - refINcigar);
		int read_end, ref_end, n_cigar_;
		uint32_t *cigar_= NULL;
		//
		//printf("ref:\t");for (i = 0; i < tlen; ++i) printf("%c", "ACGT"[(ref_seq+refINcigar)[i]]);printf("\n");
		//printf("read:\t");for (i =0; i < qlen; ++i) printf("%c", "ACGT"[(read_seq+readINcigar)[i]]);printf("\n");
		ksw_extend2(qlen, read_seq+readINcigar, tlen , ref_seq+refINcigar, 5, bwasw_sc_mat, 5, 2, abs(read_len-readINcigar-ref_len + refINcigar)+3, hash_len*bwasw_sc_mat[0], &read_end, &ref_end, &n_cigar_, &cigar_);
		if (cigar_ != NULL)
		{
			_push_cigar(cigar, cigar_len, cigar_m, cigar_, n_cigar_);
			//for (i = 0; i < n_cigar_; ++i) (*cigar)[(*cigar_len)++] = cigar_[i];
			//printf("tail:\t"); printcigar(cigar_, n_cigar_); printf("\n");
			free(cigar_);
		}
		if (read_end < read_len - readINcigar) 
			(*cigar)[(*cigar_len)++] = (((read_len - readINcigar-read_end) << 4) | CSOFT_CLIP); //'S' exist
		//(*cigar)[(*cigar_len)++] = ((read_len - readINcigar) << 4) | CSOFT_CLIP;
	}

	//printf("right bound:\t"); printcigar(*cigar, *cigar_len); printf("\n");
	free(hash_pos);
	return res;
}

int hash_left_bound_map(uint32_t **cigar, int *cigar_len, int *cigar_m,
						uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
						uint32_t **hash_num, uint64_t ***hash_node, 
						int hash_len, int hash_key, int hash_step)
{
	int hash_cigar_len=0, hash_cigar_m=100;
	uint32_t *hash_cigar = (uint32_t*)malloc(hash_cigar_m * sizeof(uint32_t));;
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	(*cigar_len) = 0;
	int hash_size = (int)pow((double)NT_N, (double)hash_key);	
	if (init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_pos, hash_key, hash_size) != 0)
	{
		(*cigar_len = 0);
		free(hash_pos);
		return 1;
	}
	int res = hash_split_map(&hash_cigar, &hash_cigar_len, &hash_cigar_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, hash_key, *hash_num, *hash_node, hash_pos, 0, 1);
	//res = hash_split_map(cigar, cigar_len, cigar_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, hash_key, *hash_num, *hash_node, hash_pos, 0, 1);
	int i, readINcigar=0, refINcigar=0;
	for (i = 0; i < hash_cigar_len; ++i)
	{
		if ((hash_cigar[i] & 0xf) == CMATCH || (*cigar)[i] & 0xf == CSOFT_CLIP) readINcigar+=(hash_cigar[i] >> 4), refINcigar+=(hash_cigar[i] >> 4);
		else if ((hash_cigar[i] & 0xf) == CINS) readINcigar += (hash_cigar[i] >> 4);
		else if ((hash_cigar[i] & 0xf) == CDEL) refINcigar += (hash_cigar[i] >> 4);
		else fprintf(stderr, "[hash_left_bound_map] cigar error: "), printcigar(*cigar, *cigar_len), fprintf(stderr, "\n"), exit(-1);
	}
	if (readINcigar < read_len)     //'S' exists
	{
		int qlen = read_len - readINcigar;
		int tlen = read_len - readINcigar + ((hash_len+read_len-readINcigar)*bwasw_sc_mat[0]-5-1)/2;
		tlen = tlen < (ref_len -refINcigar) ? tlen : (ref_len - refINcigar);
		uint8_t *qseq = (uint8_t*)malloc(qlen * sizeof(uint8_t));
		uint8_t *tseq = (uint8_t*)malloc(tlen * sizeof(uint8_t));
		int read_end, ref_end, n_cigar_;
		uint32_t *cigar_= 0;
		for (i = 0; i < qlen; ++i) qseq[i] = read_seq[qlen-i-1];
		for (i = 0; i < tlen; ++i) tseq[i] = ref_seq[ref_len-refINcigar-i-1];

		//printf("ref:\t");for (i = 0; i < tlen; ++i) printf("%c", "ACGT"[tseq[i]]);printf("\n");
		//printf("read:\t");for (i =0; i < qlen; ++i) printf("%c", "ACGT"[qseq[i]]);printf("\n");
		ksw_extend2(qlen, qseq, tlen, tseq, 5, bwasw_sc_mat, 5, 2, abs(qlen-tlen)+3, hash_len*bwasw_sc_mat[0], &read_end, &ref_end, &n_cigar_, &cigar_);

		if (cigar_ != NULL)
		{
			//printf("head:\t"); printcigar(cigar_, n_cigar_); printf("\n");
			if (read_end < read_len - readINcigar) 
				cigar_[n_cigar_++] = (((read_len -readINcigar - read_end) << 4) | CSOFT_CLIP); //'S' exsit

			//reverse_cigar
			uint32_t tmp;
			for (i = 0; i < n_cigar_/2; ++i)
				tmp = cigar_[i], cigar_[i] = cigar_[n_cigar_-i-1], cigar_[n_cigar_-i-1] = tmp;
			_push_cigar(cigar, cigar_len, cigar_m, cigar_, n_cigar_);
			free(cigar_);
		}
		else
			(*cigar)[(*cigar_len)++] = (((read_len - readINcigar) << 4) | CSOFT_CLIP); //'S' exsit

		free(tseq); free(qseq); 
	}
	_push_cigar(cigar, cigar_len, cigar_m, hash_cigar, hash_cigar_len);

	//printf("left bound:\t"); printcigar(*cigar, *cigar_len); printf("\n");
	free(hash_pos); free(hash_cigar);
	return res;
}

//for 'MIS-MATCH' case
int hash_both_bound_map(uint32_t **cigar, int *cigar_len, int *cigar_m,
						uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, 
						uint32_t **hash_num, uint64_t ***hash_node, 
						int hash_len, int hash_key, int hash_step)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	int hash_size = (int)pow((double)NT_N, (double)hash_key);	
	if (init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_pos, hash_key, hash_size) != 0)
	{
		(*cigar_len = 0);
		free(hash_pos);
		return 1;
	}
	int res;
	res = hash_split_map(cigar, cigar_len, cigar_m, ref_seq, ref_len, read_seq, read_len, hash_len, hash_step, hash_key, *hash_num, *hash_node, hash_pos, 1, 1);
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
	
	//printf("both_bound:\t"); printcigar(*cigar, *cigar_len); printf("\n");
	free(hash_pos);
	return res;
}
