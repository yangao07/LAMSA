#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "split_mapping.h"
#include "bntseq.h"
#include "lsat_aln.h"
#include "frag_check.h"
#include "ksw.h"

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
	for (i = 0; i < key_len; ++i) (*key_int) = (*key_int) << 2 | (int)seed[i];
	for (; i < hash_len; ++i) (*kmer_int) = (*kmer_int) << 2 | (int)seed[i];	
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
int init_hash_core(uint8_t *hash_seq, int hash_offset, int hash_len, int key_len, int hash_size, uint32_t *hash_num, uint64_t ***hash_node, int32_t **hash_pos)
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
//key_len = 2
//hash_size = 4^2 = 16
int init_hash(uint8_t *read_seq, int read_len, int hash_len, uint32_t **hash_num, uint64_t ***hash_node, int32_t **hash_pos, int key_len, int hash_size)
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

int hash_exact_map(uint64_t **hash_node, int key_int, int node_i, int32_t *hash_pos, int offset)
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
}

int check_hash(int h_offset, int h_i)
{
    int dis = h_offset - h_i;
    if (dis <= 10 && dis >= -10)
        return MATCH;
    else if (dis > 10)
        return DELETION;
    else return INSERT;
}

int hash_map(uint32_t *hash_num, uint64_t **hash_node, int32_t *hash_pos, int *r_i, uint8_t *query, int offset, int hash_len, int key_len)
{
	int q_key_int, q_kmer_int, node_i;

	hash_calcu(&q_key_int, &q_kmer_int, query, hash_len, key_len);
    //printf("%d %d ", q_key_int, q_kmer_int);
	//hash-hit is different with hash-search?
	if (hash_hit(hash_num, hash_node, &node_i, q_key_int, q_kmer_int, hash_len-key_len) == 1)//hit
	{
		*r_i = hash_exact_map(hash_node, q_key_int, node_i, hash_pos, offset);
		//printf("ref %d read %d \n", offset, *r_i);
        if (check_hash(offset, *r_i) == MATCH)
            return 1;
        else return 2;
	}
	else return 0;
}

int hash_get_dis(int *hash_pos, int *start_a, int a, int a_i, int b, int b_i, int hash_len, int *con_flag)
{
	if (a == b && a_i == b_i) {*con_flag = MATCH; return 0;}
	int exp, act, dis;

	/*if (abs(a-b) < 10)
	{
		if ((hash_pos[start_a[a]+a_i] - hash_pos[start_a[b]+b_i]) == (a-b))
		{
			*con_flag = MATCH; return 0;
		}
		else 
		{
			*con_flag = UNCONNECT; return -1;
		}
	}*/
	//abs(a-b) >= 10:
	exp = hash_pos[start_a[a]+a_i] + b - a;
	act = hash_pos[start_a[b]+b_i];
	dis = (abs(b-a)/(b-a)) * (act - exp);


	if (dis >= 5 && dis < DEL_THD) *con_flag = INSERT;
    else if (dis <= -5 && dis >= (0-(abs(a-b)-hash_len))) *con_flag = DELETION;
    else if (dis < 5 && dis > -5) *con_flag = MATCH; 
	else *con_flag = UNCONNECT;
	//printf("%d", dis);
	return abs(dis);
}

/*int hash_main_dis(int *hash_pos, int *start_a, int a, int a_i, int b, int b_i, int hash_len, int *con_flag)
{
	if (a == b && a_i == b_i) {*con_flag = MATCH; return 0;}
	int exp, act, dis;

	exp = hash_pos[start_a[a]+a_i] + b - a;
	act = hash_pos[start_a[b]+b_i];
	dis = (act - exp);

	if (dis == 0) *con_flag = MATCH; 
	else *con_flag = UNCONNECT;
	//printf("%d", dis);
	return abs(dis);
}*/ 
int hash_main_dis(int *hash_pos, int *start_a, int a, int a_i, int offset, int hash_len, int *con_flag)
{
	if ((hash_pos[start_a[a]+a_i] - a) == offset)
	{
		*con_flag = MATCH; 
		return 0;
	}
	else
	{
		*con_flag = UNCONNECT;
		return abs(hash_pos[start_a[a]+a_i] - a - offset);
	}
}

int hash_new_line(int offset, int hash_i, int **line, int *line_end, int path_n)
{
	////line[i][0] : the first unique-aln seed.
	//line[i][0]: the offset (read_pos - ref_pos)
	line[path_n][0] = offset;
	line[path_n][1] = hash_i;
	line_end[path_n] = 1;
	return 0;
}

int hash_copy_line(int **line, int *line_end, int from, int m_i)
{
	/*for (i = 1; i <= line_end[from]; ++i)
	{
		line[to][line_end[to]+i-1] = line[from][i];
	}
	line_end[to] += line_end[from];*/

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
}

int hash_add_line(int **line, int *line_end, int line_i, int seed_i)
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
}

hash_line_t *hash_init_line(int seq_len)
{
	hash_line_t *line;
	line = (hash_line_t*)malloc(sizeof(hash_line_t));
	line->n_blanks = 1;	//XXX initial value=1, MAX value=10.
	line->blank = (hash_blank_t*)malloc(10 * sizeof(hash_blank_t));
	line->blank[0].offset = 0;
	line->blank[0].len = seq_len;

	return line;
}

int hash_modify_blank(hash_line_t *h_line, int blank_num, int b_start, int b_end)
{
	if (blank_num >= 10)	//XXX 10
	{
		fprintf(stderr, "[hash_modify_blank] ERROR: blank number is over the range.\n"); exit(-1);
	}
	h_line->blank[blank_num].offset = b_start;
	h_line->blank[blank_num].len = b_end-b_start+1;
	return 0;
}
int hash_delete_blank(hash_line_t *h_line, int blank_num)
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
}
int hash_insert_blank(hash_line_t *h_line, int blank_num, int b_start, int b_end)
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
}
				/*ref_line | read_line*/ /*ref_len, read_len*/
int hash_fill_blank(hash_line_t *h_line, int seq_len, int blank_num, int line_start, int line_len)
{
	if (blank_num == -1) return 0;
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
}

//blank_cover: return the covered bps, set blank_i as the index of blank
//if blank_i equal -1, that indicates that this node-line doesn't cover any blank
int blank_cover(hash_line_t *h_line, int *blank_i, int pos_start, int len)
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
}

int hash_line_bound(int **line, int *line_end, int *start_a, int *hash_pos, int line_i, int hash_len, int *ref_start, int *ref_end, int *read_start, int *read_end)
{
	int offset = line[line_i][0];

	*ref_start = line[line_i][1];
	*ref_end = line[line_i][line_end[line_i]] + hash_len - 1;
	*read_start = (*ref_start) + offset;
	*read_end = (*ref_end) + offset;

	return 0;
}

//pick_the_line: find out the node-line which cover the most bps in blank, set the ref_blank_i and read_blank_i, return the line index.
//if the return value equals -1, it indicates that NO node-line covers any blank.
int pick_the_line(int **line, int *line_end, int *copy_flag, int *start_a, int *hash_pos, int hash_len, int path_n, hash_line_t *ref_line, hash_line_t *read_line, int *ref_blank_i, int *read_blank_i)
{
	int i, max_i, max_cover, cover, tmp_ref_blank_i, tmp_read_blank_i;
	int ref_start, ref_end, read_start, read_end;
	//max_cover = (int*)calloc(path_n, sizeof(int));

	max_i = -1;
	max_cover = 0;
	for (i = 0; i < path_n; ++i)
	{
		if (copy_flag[i] != 1)
		{
			hash_line_bound(line, line_end, start_a, hash_pos, i, hash_len, &ref_start, &ref_end, &read_start, &read_end);
			//line[i][1], line[i][line_end[i]] => ref_start, ref_len;
			//line[i][0], hash_pos[], start_a[] => read_start, read_len

			//tmp_blank_i == -1  means that the line can't fill any blanks.
			cover = blank_cover(ref_line, &tmp_ref_blank_i, ref_start, ref_end-ref_start+1)
				  //+ blank_cover(read_line, &tmp_read_blank_i, read_start, read_end-read_start+1);
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
}

//cut long node-line that has big gap
int hash_move_line(int **line, int *line_end, int from, int start, int end, int to)
{
	int i;

	line[to][0] = line[from][0];
	for (i = start; i <= end; ++i)
	{
		line[to][i-start+1] = line[from][i];
	}
	line_end[to] = end - start + 1;
	return 0;
}

int hash_cut_line(int **line, int *line_end, int *line_n, int hash_len)
{
	int i, j, tmp_n = *line_n;
	for (i = 0; i < tmp_n; ++i)
	{
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
}

//line: i
int hash_main_line(int *hash_pos, int *start_a, int *len_a, int ref_len, int read_len, int hash_len, int **line, int *line_end)
{
	int i, j, k, path_n = 0;
	int flag, con_flag;
	int path_i[500], tmp, copy_flag[500];

	//filter out unique-aln seeds
		for (i = 0; i <= ref_len-hash_len; ++i)
		{
			if (len_a[i] == 1)
			{
				flag = 0;
				for (j = path_n-1; j >= 0; --j)
				{
					hash_main_dis(hash_pos, start_a, i, 0, line[j][0], hash_len, &con_flag);
					if (con_flag == MATCH)
					{
						flag = 1;
						line[j][line_end[j]+1] = i;
						++line_end[j];
						if (j != (path_n-1))
						{
							for (k = j+1; k != path_n; ++k)
								line_end[k] = 0;
							path_n = j+1;
						}
						break;
					}
					//XXX
					if (line_end[j] > 5) break;
				}
				if (flag == 0)
				{
					hash_new_line(hash_pos[start_a[i]]-i, i, line, line_end, path_n);
					++path_n;
				}
			}
		}
		if (path_n < 1)
		{
			line_end[0] = 0;
			return 0;
		}

	//add multi-aln seed to the unique-aln seeds lines.
		for (i = 0; i <= ref_len-hash_len; ++i)
		{
			if (len_a[i] > 1)
			{
				for (k = 0; k < path_n; ++k)
				{
					for (j = 0; j < len_a[i]; ++j)
					{
						hash_main_dis(hash_pos, start_a, i, j, line[k][0], hash_len, &con_flag);
						if (con_flag == MATCH)
						{
							hash_add_line(line, line_end, k, i);
							break;
						}
					}
				}
				/*
				for (j = 0; j < len_a[i]; ++j)
				{
					for (k = 0; k < path_n; ++k)
					{//int hash_get_dis(int *hash_pos, int *start_a, int a, int a_i, int b, int b_i, int hash_len, int *con_flag)
						//if (i < line[line_end[path_i[k]]-1])
						hash_main_dis(hash_pos, start_a, line[path_i[k]][0], i, j, hash_len, &con_flag);
						if (con_flag == MATCH)
						{
							hash_add_line(line, line_end, path_i[k], i);
							break;
						}
					}
				}*/
			}
		}
		hash_cut_line(line, line_end, &path_n, hash_len);
		/*for (i = 0; i < path_n; ++i) 
		{
			printf("	%d: %d\n", i, line[i][0]);
			printf("%d -> %d\n", line[i][1], line[i][line_end[i]]);
		}*/

        if (line_end[0] == 59)
            printf("debug");
		fprintf(stdout, "all line:\nstart: ");
		for (i = 0; i < path_n; ++i)
		{
			printf("%d, ", line[i][1]);
		}
		printf("\nend: ");

		for (i = 0; i < path_n; ++i)
		{
			printf("%d, ", line[i][line_end[i]]+hash_len-1);
		}
		printf("\noffset: ");
		for (i = 0; i < path_n; ++i)
		{
			printf("%d, ", line[i][0]);
		}
		printf("\n");
	

    //sort the length of every line
    //XXX based on the number of nodes.
		for (i = 0; i < path_n; ++i) { path_i[i] = i; copy_flag[i] = 0; }
		for (i = 0; i < path_n-1; ++i)
		{
			for (j = i+1; j < path_n; ++j)
			{
				if(line_end[path_i[i]] < line_end[path_i[j]])
				{
					tmp = path_i[i];
					path_i[i] = path_i[j];
					path_i[j] = tmp;
				}
			}
		}
    /*
		//combine other small lines to "max_i" line
		//XXX small lines do NOT match with each other:
		//remain the long one.
		line_end[path_n] = 0;
		copy_flag[path_i[0]] = 1;

		//path_i[0] : the longest path
		//
		for (i = 1; i < path_n; ++i)
		{
			flag = 1;
			for (j = 0; j < i; ++j)
			{

				//printf("dis : %d %d -> ", path_i[i], path_i[j]);
				if (path_i[i] < path_i[j])
					hash_get_dis(hash_pos, start_a, line[path_i[j]][0], 0, line[path_i[i]][line_end[path_i[i]]-1], 1, hash_len, &con_flag);
				else
					hash_get_dis(hash_pos, start_a, line[path_i[i]][0], 0, line[path_i[j]][line_end[path_i[j]]-1], 1, hash_len, &con_flag);
				if (con_flag == UNCONNECT)
				{
					flag = 0;
					break;
				}
			}
			copy_flag[path_i[i]] = flag;
		}

		for (i = 0; i < path_n; ++i)
		{
			if (line_end[i] == 1) continue;	//delete all the line contain only one seed.
			if (copy_flag[i] == 1) hash_copy_line(line, line_end, i, path_n);
	}*/
    //greed algorithm
    int ref_start, ref_end, read_start, read_end;
    hash_line_t *ref_line, *read_line;
    ref_line = hash_init_line(ref_len); read_line = hash_init_line(read_len);
		//initialization of blanks for first line
    /*ref_start = line[path_i[0]][1];
    ref_end = line[path_i[0]][line_end[path_i[0]]]+hash_len-1;

    line_offset = hash_pos[start_a[line[path_i[0]][0]]] - line[path_i[0]][0];

    read_start = ref_start + line_offset;
    read_end = ref_end + line_offset; */
	hash_line_bound(line, line_end, start_a, hash_pos, path_i[0], hash_len, &ref_start, &ref_end, &read_start, &read_end);

    hash_fill_blank(ref_line, ref_len, 0, ref_start, ref_end-ref_start+1);
	/*fprintf(stdout, "ref:\n");
	for (i = 0; i < ref_line->n_blanks; ++i)
		fprintf(stdout, "%d -> %d\n", ref_line->blank[i].offset, ref_line->blank[i].offset+ref_line->blank[i].len-1);*/

    hash_fill_blank(read_line, read_len, 0, read_start, read_end-read_start+1);
	/*fprintf(stdout, "read:\n");
	for (i = 0; i < read_line->n_blanks; ++i)
		fprintf(stdout, "%d -> %d\n", read_line->blank[i].offset, read_line->blank[i].offset+read_line->blank[i].len-1); */
	copy_flag[path_i[0]] = 1;
	//calculate the quantity of remaining blank.
	int n_ref_blank = ref_line->n_blanks;
	int n_read_blank = read_line->n_blanks;
	int rem_line_num = path_n - 1;

	int line_i, ref_blank_i, read_blank_i;
	//XXX check nodes in lines first
	//then, if blank still exist, check nodes that are NOT in lines.
	//remaining node-lines.
    while (n_ref_blank >= 1 && n_read_blank >= 1 && rem_line_num >= 1)
    {
		//find the node-line which can cover the most of the remaining blank.
		//(line_i, blank_i) = pick_the_line();
		line_i = pick_the_line(line, line_end, copy_flag, start_a, hash_pos, hash_len, path_n, ref_line, read_line, &ref_blank_i, &read_blank_i);
		if (line_i == -1) break;
		copy_flag[line_i] = 1;
		--rem_line_num;
		//get ref_start, ref_len
		//get read_start, read_len
		hash_line_bound(line, line_end, start_a, hash_pos, line_i, hash_len, &ref_start, &ref_end, &read_start, &read_end);
		hash_fill_blank(ref_line, ref_len, ref_blank_i, ref_start, ref_end-ref_start+1);
		hash_fill_blank(read_line, read_len, read_blank_i, read_start, read_end-read_start+1);
	/*	fprintf(stdout, "ref:\n");
		for (i = 0; i < ref_line->n_blanks; ++i)
			fprintf(stdout, "%d -> %d\n", ref_line->blank[i].offset, ref_line->blank[i].offset+ref_line->blank[i].len-1);

		fprintf(stdout, "read:\n");
		for (i = 0; i < read_line->n_blanks; ++i)
			fprintf(stdout, "%d -> %d\n", read_line->blank[i].offset, read_line->blank[i].offset+read_line->blank[i].len-1);	*/
		n_ref_blank = ref_line->n_blanks;
		n_read_blank = read_line->n_blanks;
    } 

	fprintf(stdout, "main line:\nstart: ");
	for (i = 0; i < path_n; ++i)
	{
		if (copy_flag[i] == 1)
		{	//start
			printf("%d, ", line[i][1]);
		}
	}
	printf("\nend: ");

	for (i = 0; i < path_n; ++i)
	{
		if (copy_flag[i] == 1)
			printf("%d, ", line[i][line_end[i]]+hash_len-1);
	}
	printf("\noffset: ");
	for (i = 0; i < path_n; ++i)
	{
		if (copy_flag[i] == 1)
		{
			printf("%d, ", line[i][0]);
		}
	}
	printf("\n");
    line_end[path_n] = 0;
	for (i = 0; i < path_n; ++i)
	{
		if (copy_flag[i] == 1) hash_copy_line(line, line_end, i, path_n);// XXX len = 3 * line_num
	}
    free(ref_line); free(read_line);
	return path_n; 
}

int hash_add_path(int *hash_pos, int *start_a, int *len_a, int hash_len, path_msg **path, int *price_n, int start, int end, int rev)
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
				if (con_flag != UNCONNECT) {
                    last_dis = hash_get_dis(hash_pos, start_a, last_i, last_j, i, j, hash_len, &last_flag);
                    if (last_flag == MATCH)
                    {
                        path[i][j].from.x = last_i;
                        path[i][j].from.y = last_j;
                        path[i][j].flag = MATCH;
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
				if (tmp_flag != UNCONNECT)
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
							if (last_flag != UNCONNECT)
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
				if (con_flag != UNCONNECT) 
				{
					hash_get_dis(hash_pos, start_a, i, j, end, 0, hash_len, &end_flag);
					if (end_flag != UNCONNECT) {
                        last_dis = hash_get_dis(hash_pos, start_a, last_i, last_j, i, j, hash_len, &last_flag);
                        //XXX aln-res are both MATCH, then compare on the cigar_len and edit-dis.
                        /*if (last_flag == MATCH) {
                            path[i][j].from.x = last_i;
                            path[i][j].from.y = last_j;
                            path[i][j].flag = MATCH;
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
				if (tmp_flag != UNCONNECT)
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
							if (last_flag != UNCONNECT)
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
			path[back_i][back_j].flag = PATH_END;
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
}

int hash_frag_set(hash_frag *h_frag, int ref_bound, int read_bound, int flag, int hash_len)
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
	/*else	//HASH_FRAG_SEED
	{

	}*/
	return 0;
}

int hash_backtrack(path_msg **path, int *price_n, int *hash_pos, int *start_a, int n_hnode, int hash_len, hash_frag *h_frag)
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
    while (path[last_x][last_y].flag != PATH_END) {
    	if (path[last_x][last_y].flag != MATCH) {
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
}

int hash_frag_check(hash_frag h_frag, uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, uint32_t *cigar, int *cigar_len)
{
	int i;
	printf("ref:  ");
	for (i = h_frag.frag_num-1; i >= 0; --i)
	{
		printf("%4d -> %4d    ", h_frag.hfa[i].ref_start, h_frag.hfa[i].ref_end);
	}
	printf("\nread: ");
	for (i = h_frag.frag_num-1; i >= 0; --i)
		printf("%4d -> %4d    ", h_frag.hfa[i].read_start, h_frag.hfa[i].read_end);
	printf("\n");
	return 0;
}

//return the offset len of next line
int make_indel_cigar(int ref_left, int read_left, int ref_right, int read_right, int *clen, uint32_t **cigar)
{
	int dlen, ilen;
	dlen = ref_right - ref_left - 1;
	ilen = read_right - read_left - 1;
    if ((dlen > 0) && (ilen > 0))
    {
        fprintf(stderr, "[make_indel_cigar] Error: dlen: %d, ilen: %d.\n", dlen, ilen);
        exit(-1);
    }
    int len = dlen - ilen;
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
    return (dlen < ilen ? dlen : ilen); 
}

int _push_cigar(uint32_t **res_cigar, int *res_clen, uint32_t *cigar, int clen)
{
	if (clen == 0) return 0;
	int i, j;
	i = *res_clen;
	if (((i-1) >= 0) && (((*res_cigar)[i-1] & 0xf) == (cigar[0] & 0xf)))
	{
		(*res_cigar)[i-1] += ((cigar[0] >> 4) << 4);
		j = 1;
	}
	else j = 0;

	for (; j < clen; ++i,++j) {
		/*if (i == *fcmax) {
			(*fcmax) <<= 1;
			(*fcigar) = (uint32_t*)realloc(*fcigar, (*fcmax) * sizeof (uint32_t));
			if ((*fcigar) == NULL)	{fprintf(stderr, "[frag_check] Memory is not enougy.\n");exit(-1);}
		}*/	
		(*res_cigar)[i] = cigar[j];
	}
	*res_clen = i;
	return 0;
}

int hash_split_map(uint32_t **split_cigar, uint8_t *ref_seq, int ref_len, uint8_t *read_seq, int read_len, int hash_len, int key_len, uint32_t *hash_num, uint64_t **hash_node, int32_t *hash_pos)
{
	int i, read_i;//q_key_int, q_kmer_int, 
	uint8_t *ref_query;
	int cur_ref_bound, new_ref_bound, cur_read_bound, new_read_bound;
	int q_key_int, q_kmer_int, node_i;

	int *start_a = (int*)malloc((ref_len-hash_len+1) * sizeof(int));
	int *len_a = (int*)malloc((ref_len-hash_len+1) * sizeof(int));

	for (i = 0; i <= ref_len - hash_len; ++i)
	{
		ref_query = ref_seq + i;
		hash_calcu(&q_key_int, &q_kmer_int, ref_query, hash_len, key_len);
		if (hash_hit(hash_num, hash_node, &node_i, q_key_int, q_kmer_int, hash_len-key_len) == 1)
		{
			start_a[i] = ((hash_node[q_key_int][node_i] & 0xffff0000) >> 16);
			len_a[i] = (hash_node[q_key_int][node_i] & 0xffff);
		}
		else start_a[i] = len_a[i] = -1;
	}

	for (i = 0; i <= ref_len - hash_len; ++i) 
	{
		int j;
		fprintf(stdout, "%d %d -> ", i, len_a[i]);
		for (j = 0; j < len_a[i]; ++j)
			fprintf(stdout, "%d ", hash_pos[start_a[i]+j]);
		fprintf(stdout,"\n");
	}
	//main line
	int **line, *line_end, m_i;
	path_msg **path;
	int *price_n;

	line = (int**)malloc((ref_len-hash_len+1) * sizeof(int*));
	path = (path_msg**)malloc((ref_len-hash_len+1) * sizeof(path_msg*));
	for (i = 0; i <= ref_len-hash_len; ++i) {
		line[i] = (int*)malloc((ref_len-hash_len+1+1) * sizeof(int));
		path[i] = (path_msg*)malloc((read_len-hash_len+1) * sizeof(int));
	}
	line_end = (int*)malloc((ref_len-hash_len+1) * sizeof(int));
	price_n = (int*)malloc((ref_len-hash_len+1) * sizeof(int));

	//line[m_i][]: [(offset, start, len),(),()...]
	m_i = hash_main_line(hash_pos, start_a, len_a, ref_len, read_len, hash_len, line, line_end);
	if (line_end[m_i] == 0) {fprintf(stderr, "[hash_split_map] no hash-node is uniquely aligned to the read.\n"); return 0; }
	/*printf("main line:");
	for (i = 0; i < line_end[m_i]; ++i)
		printf("%d ", line[m_i][i]); printf("\n");*/

	//fill blank and SW, return the WHOLE cigar
	int _q_len, _t_len, _clen, _b_w, _score; 
	int split_clen=0;
    int line_off;
	uint32_t *_cigar=NULL;
	uint32_t *g_cigar;
	g_cigar = (uint32_t*)malloc(10 * sizeof(uint32_t));
	//1. fix the region between left bound and first line
	if (line[m_i][1] !=0 && (line[m_i][1] + line[m_i][0]) != 0)//blank exists
	{
		//split ref
		_t_len = line[m_i][1]; 
		//split read
		_q_len = line[m_i][1] + line[m_i][0];
		_b_w = abs(_t_len - _q_len);
		_score = ksw_global(_q_len, read_seq, _t_len, ref_seq, 5, &bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
		_push_cigar(split_cigar, &split_clen, _cigar, _clen);
	}
	else//add cigar
	{
		//         ref_left, read_left, ref_rigth,    read_right
		make_indel_cigar(-1, -1,        line[m_i][1], line[m_i][1]+line[m_i][0], &_clen, &g_cigar);
		_push_cigar(split_cigar, &split_clen, g_cigar, _clen);
	}
	//2. fix the region betweed each lines
	for (i = 0; i < line_end[m_i]; ++i)
	{
		//generate cigar for every main-lines
		//generate_mem_cigar(,);
		g_cigar[0] = ((line[m_i][i*3+2]+hash_len-1) << 4) | CMATCH; _clen = 1;
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
			_push_cigar(split_cigar, &split_clen, _cigar, _clen);
		}
		else //add cigar
		{
			line_off = make_indel_cigar((line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2), line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3], line[m_i][(i+1)*3+1], line[m_i][(i+1)*3+1]+line[m_i][(i+1)*3], &_clen, &g_cigar);
            //for overlaps between lines
            line[m_i][(i+1)*3+1] -= line_off;
            line[m_i][(i+1)*3+2] += line_off;
			_push_cigar(split_cigar, &split_clen, g_cigar, _clen);
		}
	}
	//3. fix the region between last line and right bound	//i = line_end[m_i]-1
	if ((line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2) < (ref_len - 1) && (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]) < (read_len-1))	//blank exists
	{
		_t_len = (ref_len-1) - (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2);
		_q_len = (read_len-1) - (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]);
		_score = ksw_global(_q_len, read_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len+line[m_i][i*3]-1), _t_len, ref_seq+(line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-1), 5, &bwasw_sc_mat, 5, 2, _b_w, &_clen, &_cigar);
		_push_cigar(split_cigar, &split_clen, _cigar, _clen);
	}
	else //add cigar
	{
		make_indel_cigar((line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2), (line[m_i][i*3+1]+line[m_i][i*3+2]+hash_len-2+line[m_i][i*3]), ref_len, read_len, &_clen, &g_cigar);
		_push_cigar(split_cigar, &split_clen, g_cigar, _clen);
	}
	
	if (_cigar != NULL) free(_cigar); 
	free(g_cigar);
	return split_clen;
	
	//add path
		/*if (hash_add_path(hash_pos, start_a, len_a, hash_len, path, price_n, 0, line[m_i][0], 1) == 0) //rev-path one endpoint
			path[line[m_i][0]][0].flag = PATH_END;
		for (i = 0; i < line_end[m_i]-1; i++)
		{
			//fprintf(stderr, "debug: %d %d %d %lld %lld %lld \n", i, line[m_i][i], line[m_i][i+1], a_msg[line[m_i][i]].read_id, a_msg[line[m_i][i+1]].read_id, a_msg[line[m_i][i]].at[0].offset);
			hash_add_path(hash_pos, start_a, len_a, hash_len, path, price_n, line[m_i][i], line[m_i][i+1], -2);	//two endpoint
		}
		hash_add_path(hash_pos, start_a, len_a, hash_len, path, price_n, line[m_i][line_end[m_i]-1], ref_len-hash_len, -1);	//one endpoint

		//backtrack
		hash_frag h_frag;
		h_frag.frag_num = 0;
		h_frag.frag_max = (ref_len-hash_len+1);
		h_frag.hfa = (hash_frag_aln*)malloc(h_frag.frag_max * sizeof(hash_frag_aln));

		hash_backtrack(path, price_n, hash_pos, start_a, ref_len-hash_len+1, hash_len, &h_frag);

		//generet cigar
		uint32_t *cigar;
		int cigar_len;
		hash_frag_check(h_frag, ref_seq, ref_len, read_seq, read_len, cigar, &cigar_len);

		free(start_a); free(len_a);	
		for (i = 0; i <= ref_len-hash_len; ++i) { free(line[i]); free(path[i]); }
		free(line); free(line_end); free(price_n);
		free(h_frag.hfa);*/
}


//XXX for INS: the ref and read bound maybe the right bound of INS seqs, cause of the INS is just a repeat of seqs nearby.
int find_bound(uint8_t *ref_seq, int ref_len, int read_len, int hash_len, int key_len, uint32_t *hash_num, uint64_t **hash_node, int32_t *hash_pos, int *ref_b1, int *ref_b2, int *read_b1, int *read_b2)
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
}

int separate_cigar(uint32_t *cigar, int cigar_len, int slen, uint32_t *c1, int *c1_len, uint32_t *c2, int *c2_len)
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
}

int amend_bound_core(uint32_t *c1, int c1_len, uint32_t *c2, int c2_len, int *sseq_b1, int *sseq_b2, int *lseq_b1, int *lseq_b2)
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
}

int amend_bound(uint8_t *sseq, int *sseq_b1, int *sseq_b2, uint8_t *lseq, int *lseq_b1, int *lseq_b2)
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
	if (qlen == 350 && tlen == 314)
		printf("debug");
	score = ksw_global(qlen, query, tlen, target, 5, &bwasw_sc_mat, 5, 2, b_w, &cigar_len, &cigar); printcigar(cigar, cigar_len);
	c1 = (uint32_t*)malloc(cigar_len * sizeof(uint32_t)); c2 = (uint32_t*)malloc(cigar_len * sizeof(uint32_t));
	separate_cigar(cigar, cigar_len, tlen, c1, &c1_len, c2, &c2_len);

	printcigar(c1, c1_len);	printcigar(c2, c2_len);
	amend_bound_core(c1, c1_len, c2, c2_len, sseq_b1, sseq_b2, lseq_b1, lseq_b2);

	if (qlen == 2 * tlen) free(query);
	free(cigar); free(c1); free(c2);
	return 0;
}

int get_res_cigar(uint8_t *sseq, int slen, int sseq_b1, int sseq_b2, uint8_t *lseq, int llen, int lseq_b1, int lseq_b2, uint32_t **res_cigar, int *cigar_len, int FLAG)
{
	int8_t cin;
	if (FLAG == 1) cin = CINS; else cin = CDEL;
	
	int tlen, qlen, b_w, score, c_len, c1_len, c2_len, i;
	uint32_t *c1, *c2;
	tlen = sseq_b1+1; qlen = lseq_b1+1; b_w = abs(tlen-qlen);
	if (qlen == 350 && tlen == 314)
		printf("debug");
	score = ksw_global(qlen, lseq, tlen, sseq, 5, &bwasw_sc_mat, 5, 2, b_w, &c1_len, &c1); //printcigar(c1, c1_len);
	tlen = slen - sseq_b2;	qlen = llen - lseq_b2; b_w = abs(tlen -qlen);

	if (qlen == 350 && tlen == 314)
		printf("debug");
	score = ksw_global(qlen, lseq+lseq_b2, tlen, sseq+sseq_b2, 5, &bwasw_sc_mat, 5, 2, b_w, &c2_len, &c2); //printcigar(c2, c2_len);

	//*res_cigar = (uint32_t*)malloc((c1_len+c2_len+2)*sizeof(uint32_t));
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
}

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
int split_insert_map(uint32_t **res_cigar, int *res_len, uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, int64_t ref_offset, int hash_len, uint32_t **hash_num, uint64_t ***hash_node, int key_len, int hash_size)
{
	//fprintf(stdout, "%lld ", (long long)ref_offset);
	int ref_bound1, ref_bound2, read_bound1, read_bound2;
	int i;
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));

	for (i = 0; i < read_len; ++i) printf("%d", read_seq[i]);
		printf("\n");
	for (i = 0; i < ref_len; ++i) printf("%d", ref_seq[i]);
		printf("\n");
	//init of hash table
	init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_pos, key_len, hash_size);
	
	//for multi-breakpoints
	(*res_len) = hash_split_map(res_cigar, ref_seq, ref_len, read_seq, read_len, hash_len, key_len, *hash_num, *hash_node, hash_pos);
	printcigar(*res_cigar, *res_len); fprintf(stdout,"\n");
	// bound1 and bound2
	/*find_bound(ref_seq, ref_len, read_len, hash_len, key_len, *hash_num, *hash_node, hash_pos, &ref_bound1,&ref_bound2, &read_bound1, &read_bound2);
	fprintf(stdout, "bound : %d %d %d %d ", ref_bound1, read_bound1, ref_bound2, read_bound2);
	if ((ref_bound2 - ref_bound1) > 20)
	{
		fprintf(stdout, "multi-break candi-pos.\n");
	}
	else {
		//insert allowed SW and amend bounds
		amend_bound(ref_seq, &ref_bound1, &ref_bound2, read_seq, &read_bound1, &read_bound2);
		fprintf(stdout, "amendatory bound : %d %d %d %d\n", ref_bound1, read_bound1, ref_bound2, read_bound2);
		//combine result, print split-mapping cigar
		get_res_cigar(ref_seq, ref_len, ref_bound1, ref_bound2, read_seq, read_len, read_bound1, read_bound2, res_cigar, res_len, 1);
		printcigar(*res_cigar, *res_len);
		printf("\n");
	}*/
	//printcigar(*res_cigar, *res_len);
	//fprintf(stdout, "\n\n");
	free(hash_pos);
	//exit(-1);
	return 0;
}

//deletion length = ref_len - read_len
int split_delete_map(uint32_t **res_cigar, int *res_len, uint8_t *read_seq, int read_len, uint8_t *ref_seq, int ref_len, int64_t ref_offset, int hash_len, uint32_t **hash_num, uint64_t ***hash_node, int key_len, int hash_size)
{
	int32_t *hash_pos = (int32_t*)malloc(read_len * sizeof(int32_t));
	//fprintf(stdout, "%lld ", (long long)ref_offset);
	int ref_bound1, ref_bound2, read_bound1, read_bound2;
	int i;
	for (i = 0; i < read_len; ++i) printf("%d", read_seq[i]);
		printf("\n");
	for (i = 0; i < ref_len; ++i) printf("%d", ref_seq[i]);
		printf("\n");
	//init of hash table
	init_hash(read_seq, read_len, hash_len, hash_num, hash_node, &hash_pos, key_len, hash_size);

	//for multi-breakpoints
	(*res_len) = hash_split_map(res_cigar, ref_seq, ref_len, read_seq, read_len, hash_len, key_len, *hash_num, *hash_node, hash_pos);
	printcigar(*res_cigar, *res_len); fprintf(stdout,"\n");

	// bound1 and bound2
	/*find_bound(ref_seq, ref_len, read_len, hash_len, key_len, *hash_num, *hash_node, hash_pos, &ref_bound1, &ref_bound2, &read_bound1, &read_bound2);
	fprintf(stdout, "bound : %d %d %d %d ", ref_bound1, read_bound1, ref_bound2, read_bound2);

	if ((read_bound2 - read_bound1) > 20)
	{
		fprintf(stdout, "multi-break candi-pos.\n");
	}
	else {
		//deletion allowed SW XXX+
		amend_bound(read_seq, &read_bound1, &read_bound2, ref_seq, &ref_bound1, &ref_bound2);
		fprintf(stdout, "amendatory bound : %d %d %d %d\n", ref_bound1, read_bound1, ref_bound2, read_bound2);

		//combine result, print split-mapping cigar

		get_res_cigar(read_seq, read_len, read_bound1, read_bound2, ref_seq, ref_len, ref_bound1, ref_bound2, res_cigar, res_len, 0);
		printcigar(*res_cigar, *res_len);
		printf("\n");
	}*/
	//printcigar(*res_cigar, *res_len);
	//fprintf(stdout, "\n\n");
	free(hash_pos);
	//exit(-1);
	return 0;
}
