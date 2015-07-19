#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lsat_aln.h"
#include "split_mapping.h"
#include "lsat_heap.h"

//XXX 在输出node时就已经排序？即在dp过程中排序？
void line_sort_endpos(line_node **line, int *line_end,
                      trig_node **t_node, int *tri_n,
                      int li, int len) {
    int i, j, k, end, t;
    line_node tmp; trig_node tr;
    for (i = li; i < li+len-1; ++i) {
        for (j = i+1; j < len; ++j) {
            if (line[j][line_end[j]-1].x > line[i][line_end[i]-1].x) {
                // swap_line
                end = line_end[j] > line_end[i] ? line_end[j] : line_end[i];
                for (k = 0; k < end+L_EXTRA; ++k) {
                    tmp = line[i][k], line[i][k] = line[j][k], line[j][k] = tmp;   
                }
                t = line_end[i], line_end[i] = line_end[j], line_end[j] = t;
                end = tri_n[j] > tri_n[i] ? tri_n[j] : tri_n[i];
                for (k = 0; k < end; ++k) {
                    tr = t_node[i][k], t_node[i][k] = t_node[j][k], t_node[j][k] = tr;   
                }
                t = tri_n[i], tri_n[i] = tri_n[j], tri_n[j] = t;
            }
        }
    }
}
void line_sort_endpos1(line_node **line, int *line_end,
                       int li, int len) {
    int i, j, k, end, t;
    line_node tmp;
    for (i = li; i < li+len-1; ++i) {
        for (j = i+1; j < len; ++j) {
            if (line[j][line_end[j]-1].x > line[i][line_end[i]-1].x) {
                // swap_line
                end = line_end[j] > line_end[i] ? line_end[j] : line_end[i];
                for (k = 0; k < end+L_EXTRA; ++k) {
                    tmp = line[i][k], line[i][k] = line[j][k], line[j][k] = tmp;   
                }
                t = line_end[i], line_end[i] = line_end[j], line_end[j] = t;
            }
        }
    }
}

node_score *node_init_score(int n) {
    node_score *ns = (node_score*)malloc(sizeof(node_score));
    ns->max_n = n;
    ns->node_n = 0;
    ns->node = (line_node*)malloc(n * sizeof(line_node));
    ns->score = (int*)malloc(n * sizeof(int));
    return ns;
}

void node_free_score(node_score *ns) {
    free(ns->score); free(ns->node);
    free(ns);
}

int heap_add_node(node_score *ns, line_node node, int score)
{
    if (ns->node_n < ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->node[ns->node_n++] = node;
        return -1;
    } else if (ns->node_n == ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->node[ns->node_n++] = node;
        build_node_min_heap(ns);
        return -1;
    } else { // (ns->node_n == max_n, min-heap is already built.
        return node_heap_update_min(ns, node, score); 
    }
}

int line_merge(int a, int b, line_node **line, int *line_end) {
    int s1, e1, s2, e2, s, e, hi;
    float rat1, rat2;

    s2 = line[a][0].x, e2 = line[a][line_end[a]-1].x;
    //if (line[b][line_end[b]+1].x == -1) { // NO merged-line
    if (L_MF(line,line_end,b) & L_NMERG) {// Not merged
        hi = b;
        s1 = line[b][0].x, e1 = line[b][line_end[b]-1].x;
    } else { // [b] was merged already
        hi = L_MH(line, line_end, b);
        s1 = L_LB(line, line_end, hi), e1 = L_RB(line, line_end, hi);
    }
    s = (s2 > s1) ? s2 : s1;
    e = (e2 < e1) ? e2 : e1;
	/*
    // new strategy: no mattter how big/small the overlap is, do merge.
    if (s <= e) {
        L_LB(line, line_end, hi) = s1+s2-s, L_RB(line, line_end, hi) = e1+e2-e;
        L_MF(line, line_end, hi) = L_MERGH;
        L_MF(line, line_end, a) = L_MERGB, L_MH(line, line_end, a) = hi;
        return 1;
    } else {
        L_MF(line, line_end, a) = L_NMERG;
        return 0;
    }*/
    
    rat1 = (e-s+1+0.0)/(e1-s1+1+0.0);
    rat2 = (e-s+1+0.0)/(e2-s2+1+0.0);
    //if ((rat1>=0.7 || rat2>=0.7)) { //0.7 XXX
	if (rat1 < 0.7 && rat2 < 0.7) {
			L_MF(line, line_end, a) = L_NMERG;
			return 0;
	} else if 
		(L_LS(line, line_end, a) <= L_LS(line, line_end, b) / 2 || // merged and  dumped
		 L_LS(line, line_end, a) <= L_BS(line, line_end, b) / 2 ) {

        L_LB(line, line_end, hi) = s1, L_RB(line, line_end, hi) = e1;
        L_MF(line, line_end, hi) = L_MERGH;
        L_MF(line, line_end, a) = L_MERGB, L_MH(line, line_end, a) = hi;
		L_MF(line, line_end, a) |= L_DUMP;
		return 1;
	} else { // merged
        //line[hi][line_end[hi]+1].x = s1+s2-s, line[hi][line_end[hi]+1].y = e1+e2-e;
        L_LB(line, line_end, hi) = s1+s2-s, L_RB(line, line_end, hi) = e1+e2-e;
        L_MF(line, line_end, hi) = L_MERGH;
        //line[a][line_end[a]+1].x = -2, line[a][line_end[a]+1].y = hi;
        L_MF(line, line_end, a) = L_MERGB, L_MH(line, line_end, a) = hi;
		if (L_BS(line, line_end, b) > L_BS(line, line_end, a))
			L_BS(line, line_end, a) =  L_BS(line, line_end, b);
        return 1;
    } 
}

// select best/secondary merged-line, based on line-node number OR based on line-score? XXX
// rule: 1. Better Score
//       2. Short line, if Score is equal
//            
// Best & Secondary: 1. Best is unique, 
//                   2. Secondary is allowed to be multi, 
//            		 3. Secondary is allowed to equal the Best
// new
// Best and all secondary whose score exceed best/2
void line_filter(line_node **line, int *line_end, int li, int len, trig_node **tri_node, int *tri_n, frag_dp_node **f_node, aln_msg *a_msg, int per_max_multi)
{
    int i, ii, j, k, l, m_ii;
    line_node **m_b = (line_node**)malloc(len * sizeof(line_node*)); // (i, score)
	int **m_f = (int**)malloc(len * sizeof(int*)); // merged line (after filtering)
	for (i = 0; i < len; ++i) {
		m_b[i] = (line_node*)malloc(len * sizeof(line_node));
		m_f[i] = (int*)malloc((len+1) * sizeof(int));
	}
	int *m_bn = (int*)malloc(len * sizeof(int)); 
	int *m_fn = (int*)malloc(len * sizeof(int));
    int m_i = -1;
    int b_score, s_score;
	int m_head, f_merge, min_l, max_r, head;
    // check merged-lines
    for (i = li; i < li+len; ++i) {
		if (L_MF(line, line_end, i) & L_DUMP) continue;
        if (L_MF(line, line_end, i) & L_NMERG) { // NOT merged
			++m_i;
			m_b[m_i][0].x = i, m_b[m_i][0].y = -2;
			m_bn[m_i] = 1;
		} else if (L_MF(line, line_end, i) & L_MERGH) { // merged, head
			++m_i;
            m_b[m_i][0].x = i, m_b[m_i][0].y = L_LS(line, line_end, i); // NO secondary
			m_bn[m_i] = 1;
		} else { // merged, body : inter-line candidate
            //if (m_i < 0) { fprintf(stderr, "[line_filter] %s BUG.\n", READ_NAME); exit(1); }
			m_b[m_i][m_bn[m_i]].x = i, m_b[m_i][m_bn[m_i]].y = L_LS(line, line_end, i);
			m_bn[m_i]++;
        }
    }
	// select best/secondary line, set flag, set inter-line
    node_score *ns = node_init_score(per_max_multi);
	for (i = 0; i <= m_i; ++i) {
		if (m_b[i][0].y == -2) { // not merged
			m_f[i][0] = m_b[i][0].x;
			m_fn[i] = 1;
			continue;
		}
		b_score = s_score = 0;
		for (j = 0; j < m_bn[i]; ++j) {
			if (m_b[i][j].y > b_score) {
				s_score = b_score;
				b_score = m_b[i][j].y;
			} else if (m_b[i][j].y > s_score)
				s_score = m_b[i][j].y;
		}
		// half ? XXX
		if (s_score > b_score/2) f_merge = 1; // get M largest scores
		else f_merge = 0;                     // get the largest score

		m_fn[i] = 1; // m_f [0: best, 1: head, 1...: body
		if (f_merge) {
            //1. init_heap
            ns->node_n = 0;
			for (j = 0; j < m_bn[i]; ++j) {
                //5. output_M_nodes
				if (m_b[i][j].y > b_score/2) {
                    //2. add_node & build_heap & update_heap
                    int ret = heap_add_node(ns, (line_node){m_b[i][j].x,-1}, m_b[i][j].y);
                    if (ret == -2) { // dumped
                        L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
                        tri_n[m_b[i][j].x] = 0;
                    } else if (ret != -1) { // ret is dumped
                        L_MF(line, line_end, ret) |= L_DUMP;
                        tri_n[ret] = 0;
                    }
                } else {
					L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
					tri_n[m_b[i][j].x] = 0;
				}
			} 	
            //3. output_M_nodes, set head/body flags
            build_node_minpos_heap(ns);
            // head
            line_node head_n = node_heap_extract_minpos(ns);
            m_head = head_n.x;
            L_MF(line, line_end, m_head) = L_MERGH;

            min_l = line[m_head][0].x;
            max_r = line[m_head][line_end[m_head]-1].x;
            if (L_LS(line, line_end, m_head) == b_score) m_f[i][0] = m_head; // best
            m_f[i][m_fn[i]++] = m_head;

            int body;
            while ((body = node_heap_extract_minpos(ns).x) != -1) {
                L_MF(line, line_end, body) = L_MERGB;
                L_MH(line, line_end, body) = m_head;
                min_l = MINOFTWO(min_l, line[body][0].x);
                max_r = MAXOFTWO(max_r, line[body][line_end[body]-1].x);
                if (L_LS(line, line_end, body) == b_score) m_f[i][0] = body; // best
                m_f[i][m_fn[i]++] = body;
            }
            // set left/right boundary
			L_LB(line, line_end, m_head) = MINOFTWO(line[m_head][0].x, min_l);
			L_RB(line, line_end, m_head) = MAXOFTWO(line[m_head][line_end[m_head]-1].x, max_r);
		} else {
			for (j = 0; j < m_bn[i]; ++j) {
				if (m_b[i][j].y == b_score) {
					L_MF(line, line_end, m_b[i][j].x) =  L_NMERG;
					m_f[i][0] = m_b[i][j].x; // best
					m_f[i][m_fn[i]++] = m_b[i][j].x;
				} else {
					L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
					tri_n[m_b[i][j].x] = 0;
				}
			}
		}
		for (ii = 1; ii < m_fn[i]; ++ii) {
			j = m_f[i][ii];
			for (k = 0; k < tri_n[j]; ++k) {
                head = -1;
				// check inter-line: from j+1 -> merge_end
                for (l = j+1; l < li+len && (L_MF(line, line_end, l) & 0x3)==0; ++l) { // L_MERGB
                    if (line[l][0].x > tri_node[j][k].n1.x && line[l][line_end[l]-1].x < tri_node[j][k].n2.x) {
                        // check for MIS(INV)
						if (f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].match_flag == F_MISMATCH ||
                                f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].match_flag == F_LONG_MISMATCH) {
                            int x_s = line[l][0].x, y_s = line[l][0].y;
                            int x_e = line[l][line_end[l]-1].x, y_e = line[l][line_end[l]-1].y;
                            int x1 = tri_node[j][k].n1.x, y1 = tri_node[j][k].n1.y, x2 = tri_node[j][k].n2.x, y2 = tri_node[j][k].n2.y;
                            int nsrand = a_msg[x_s].at[y_s].nsrand;
                            if (nsrand == a_msg[x1].at[y1].nsrand ||
                                    a_msg[x_s].at[y_s].chr != a_msg[x1].at[y1].chr ||
                                    nsrand * a_msg[x_s].at[y_s].offset < nsrand * a_msg[x2].at[y2].offset ||
                                    nsrand * a_msg[x_e].at[y_e].offset > nsrand * a_msg[x1].at[y1].offset)
                                continue;

                            E_LB(line, line_end, l) = tri_node[j][k].n1.x;
                            E_RB(line, line_end, l) = tri_node[j][k].n2.x;
                            L_MF(line, line_end, l) = L_INTER;
                            if (tri_node[j][k].tri_flag == 0x2 && (f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].trg_n & 0x2)) // pre
                                f_node[tri_node[j][k].n2.x][tri_node[j][k].n2.y].trg_n -= 0x2;
                            else if (tri_node[j][k].tri_flag == 0x1 && (f_node[tri_node[j][k].n1.x][tri_node[j][k].n1.y].trg_n & 0x1)) // next
                                f_node[tri_node[j][k].n1.x][tri_node[j][k].n1.y].trg_n -= 0x1; 
                            if (head == -1) {
                                L_MF(line, line_end, l) |= L_NMERG;
                                head = l;
                            } else {
                                L_MF(line, line_end, l) |= L_MERGB;
                                L_MH(line, line_end, l) = head;
                                L_MF(line, line_end, head) = L_INTER|L_MERGH;
                            }
                        }
                    }
                }
            }
		}
	}
    node_free_score(ns);

    // remove candidate "trigger-line", left-most or right-most line
    if (m_i > 0) {
        int L;
        // left-boundary tri-line
        m_ii = 0;
		// best line in merged-lines
        l = line[m_f[m_ii][0]][line_end[m_f[m_ii][0]]-1].x - line[m_f[m_ii][0]][0].x;
        L = line[m_f[m_ii+1][0]][line_end[m_f[m_ii+1][0]]-1].x - line[m_f[m_ii+1][0]][0].x;
        if (l < 2 && L >= 2) { // DUMP
			for (i = 1; i < m_fn[m_ii]; ++i) {
				L_MF(line, line_end, m_f[m_ii][i]) = L_DUMP;
                for (j = li; j < li+len; ++j) {
                    if (!(L_MF(line, line_end, j) & L_NMERG) && !(L_MF(line, line_end, j) & L_MERGH) && !(L_MF(line, line_end, j) & L_DUMP) && (L_MH(line, line_end, j) == m_f[m_ii][i]))
                        L_MF(line, line_end, j) = L_DUMP;
                }
            }
        }
        // right-boundary tri-line
        m_ii = m_i;
        l = line[m_f[m_ii][0]][line_end[m_f[m_ii][0]]-1].x - line[m_f[m_ii][0]][0].x;
        L = line[m_f[m_ii-1][0]][line_end[m_f[m_ii-1][0]]-1].x - line[m_f[m_ii-1][0]][0].x;
        if (l < 2 && L >= 2) { // DUMP
            for (i = 1; i < m_fn[m_ii]; ++i) {
                L_MF(line, line_end, m_f[m_ii][i]) = L_DUMP;
                for (j = li; j < li+len; ++j) {
                    if (!(L_MF(line, line_end, j) & L_NMERG) && !(L_MF(line, line_end, j) & L_MERGH) && !(L_MF(line, line_end, j) & L_DUMP) && (L_MH(line, line_end, j) == m_f[m_ii][i]))
                        L_MF(line, line_end, j) = L_DUMP;
                }
            }
        }
    }
	for (i = 0; i < len; ++i) { free(m_b[i]); free(m_f[i]); }
    free(m_b); free(m_bn); free(m_f); free(m_fn);
}

void line_filter1(line_node **line, int *line_end, int li, int len, int per_max_multi)
{
    int i,j;
	int *m_n = (int*)malloc(len * sizeof(int));
    line_node **m_b = (line_node**)malloc(len * sizeof(line_node*)); // (i, score)
	for (i = 0; i < len; ++i) m_b[i] = (line_node*)malloc(len * sizeof(line_node));
    int m_i = -1;
    int b_score, s_score;
    for (i = li; i < li+len; ++i) {
		if ((L_MF(line, line_end, i) & L_DUMP) || (L_MF(line, line_end, i) & L_NMERG)) continue; // dumped and  NOT merged
        else if (L_MF(line, line_end, i) & L_MERGH) {    // merged, head
            ++m_i;
            m_b[m_i][0].x = i, m_b[m_i][0].y = L_LS(line, line_end, i);
			m_n[m_i] = 1;
        } else { // merged, body
            //if (m_i < 0) { fprintf(stderr, "[line_filter] %s BUG.\n", READ_NAME); exit(1); }
			m_b[m_i][m_n[m_i]].x = i, m_b[m_i][m_n[m_i]].y = L_LS(line, line_end, i);
			m_n[m_i]++;
        }
    }
    int min_l, max_r, f_merge, m_head;
    node_score *ns = node_init_score(per_max_multi);
	for (i = 0; i <= m_i; ++i) { // m_b[i][0 ... m_n[i]-1]
		// get best & secondary score
		b_score = s_score = 0;
		for (j = 0; j < m_n[i]; ++j) {
			if (m_b[i][j].y > b_score) {
				s_score = b_score;
				b_score = m_b[i][j].y;
			} else if (m_b[i][j].y > s_score)
				s_score = m_b[i][j].y;
		}
		if (s_score > b_score/2) f_merge = 1;
		else f_merge = 0;
		// select best & secondary line, set flag
		if (f_merge) {
            // 1. init heap
            ns->node_n = 0;
			for (j = 0; j < m_n[i]; ++j) {
				if (m_b[i][j].y > b_score/2) {
                    // 2. add node & build heap & update heap
                    int ret = heap_add_node(ns, (line_node){m_b[i][j].x, -1}, m_b[i][j].y);
                    if (ret == -2) { // dumped
                        L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
                    } else if (ret != -1) { // ret is dumped
                        L_MF(line, line_end, ret) |= L_DUMP;
                    }
                } else L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
			} 	
            // 3. ouput_M_nodes, set head/body flags
            build_node_minpos_heap(ns);
            line_node head_n = node_heap_extract_minpos(ns);
            m_head = head_n.x;
            L_MF(line, line_end, m_head) = L_MERGH;

            min_l = line[m_head][0].x;
            max_r = line[m_head][line_end[m_head]-1].x;

            int body;
            while ((body = node_heap_extract_minpos(ns).x) != -1) {
                L_MF(line, line_end, body) = L_MERGB;
                L_MH(line, line_end, body) = m_head;
                min_l = MINOFTWO(min_l, line[body][0].x);
                max_r = MAXOFTWO(max_r, line[body][line_end[body]-1].x);
            }
            // set left/right boundary
			L_LB(line, line_end, m_head) = MINOFTWO(line[m_head][0].x, min_l);
			L_RB(line, line_end, m_head) = MAXOFTWO(line[m_head][line_end[m_head]-1].x, max_r);
		} else {
			for (j = 0; j < m_n[i]; ++j) {
				if (m_b[i][j].y == b_score)
					L_MF(line, line_end, m_b[i][j].x) =  L_NMERG;
				else L_MF(line, line_end, m_b[i][j].x) |= L_DUMP;
			}
		}
	}
    node_free_score(ns);
	for (i = 0; i < len; ++i) free(m_b[i]);
	free(m_b); free(m_n);
}

int line_remove(line_node **line, int *line_end, int li, int len)
{
    int *head = (int*)malloc(len * sizeof(int));
    int l, j, cur_i=li;
    // remove line with "DUMP" flag
    for (l = li; l < li+len; ++l) {
        if (!(L_MF(line, line_end, l) & L_DUMP)) {
            if (L_MF(line, line_end, l) & L_MERGH)
                head[l-li] = cur_i;
            if (l == cur_i) {
                ++cur_i;
                continue;
            }
            line_end[cur_i] = line_end[l];
            for (j = 0; j < line_end[l]+L_EXTRA; ++j)
                line[cur_i][j] = line[l][j];
            //merged, body
            if (!(L_MF(line, line_end, cur_i)&L_NMERG) && !(L_MF(line, line_end, cur_i)&L_MERGH))
                    L_MH(line, line_end, cur_i) = head[L_MH(line, line_end, cur_i)-li];
            ++cur_i;
        }
    }
    free(head);
    return cur_i-li;
}


//XXX 将NotMerged与Merged统一起来
/* line: [0 ~ end-1] -> nodea; 
         [end] -> (left-b, right-b); 
         [end+1] -> (merged-flag, merge-head)/(start, end) 
         // (-1,):NOT merged
         // (-2,i):Merged and head is i
         // (-2,t):inter-merged, head is t
         // (-3,x):Merged, but dumped
         // (a,b): head, start from a to b
         [end+2] -> (line-score, x)
         */
int line_set_bound(line_node **line, int *line_end, 
                   int li, int *o_len, int left, int right, 
                   trig_node **t_node, int *tri_n, frag_dp_node **f_node, aln_msg *a_msg, int per_max_multi) 
{
    if (*o_len <= 0) return 0;
    int i, j, len = *o_len;
    // sort by end-pos
    line_sort_endpos(line, line_end, t_node, tri_n, li, len);
    // merge lines
    ///line[li][line_end[li]+1].x = -1;
    L_MF(line, line_end, li) = L_NMERG;
    for (i = 1; i < len; ++i)
        line_merge(li+i, li+i-1, line, line_end);
    // filter best/secondary line, and their tri-nodes, then find inter-line 
    line_filter(line, line_end, li, len, t_node, tri_n, f_node, a_msg, per_max_multi);
    *o_len = line_remove(line, line_end, li, len);
    len = *o_len;
    // set line-boundary
    int l,r,m,s,e;
    r = right;
    for (i = li; i < li+len; ++i) {
        if (!(L_MF(line, line_end, i) & L_DUMP)) {
            E_RB(line, line_end, i) = r;
            break;
        }
    }
    //if (i == li+len) {fprintf(stderr, "[line_set_bound] %s BUG 0.\n", READ_NAME); exit(1);}
    // set right-b
    if (L_MF(line, line_end, i) & L_NMERG)
        s = line[i][0].x, e = line[i][line_end[i]-1].x;
    else
        s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
    m = 0; //Merged head OR NOT-merged head
    for (++i; i < li+len; ++i) {
        if (L_MF(line, line_end, i) & (L_DUMP|L_INTER)) continue;
        if (L_MF(line, line_end, i) & L_NMERG) {
            r = s, l = line[i][line_end[i]-1].x;
            //set left-b for last-line(s)
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = line[i][0].x, e = line[i][line_end[i]-1].x;
            m = i;
        } else if (L_MF(line, line_end, i) & L_MERGH) {
            r = s, l = L_RB(line, line_end, i);
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
            m = i;
        } //merged and body
        E_RB(line, line_end, i) = r;
    }
    for (j = i-1; j >= m; --j) {
        if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
        E_LB(line, line_end, j) = left;
    }
    // set inter-bound
    //line_inter_bound(line, line_end, li, len, t_node, tri_n, f_node, a_msg);
    return 0;
}

int line_set_bound1(line_node **line, int *line_end, 
                   int li, int *o_len, int left, int right, int per_max_multi) 
{
    if (*o_len <= 0) return 0;
    int i, j, len=*o_len;
    // sort by end-pos
    line_sort_endpos1(line, line_end, li, len);
    // merge lines
    ///line[li][line_end[li]+1].x = -1;
    L_MF(line, line_end, li) = L_NMERG;
    for (i = 1; i < len; ++i)
        line_merge(li+i, li+i-1, line, line_end);
    // filter best/secondary line, and their tri-nodes, then find inter-line 
    line_filter1(line, line_end, li, len, per_max_multi);
    *o_len = line_remove(line, line_end, li, len);
    len = *o_len;
    // set line-boundary
    int l,r,m,s,e;
    r = right;
    for (i = li; i < li+len; ++i) {
        if (!(L_MF(line, line_end, i) & L_DUMP)) {
            E_RB(line, line_end, i) = r;
            break;
        }
    }
    //if (i == li+len) {fprintf(stderr, "[line_set_bound] %s BUG 0.\n", READ_NAME); exit(1);}
    // set right-b
    if (L_MF(line, line_end, i) & L_NMERG)
        s = line[i][0].x, e = line[i][line_end[i]-1].x;
    else
        s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
    m = 0; //Merged head OR NOT-merged head
    for (++i; i < li+len; ++i) {
        if (L_MF(line, line_end, i) & (L_DUMP|L_INTER)) continue;
        if (L_MF(line, line_end, i) & L_NMERG) {
            r = s, l = line[i][line_end[i]-1].x;
            //set left-b for last-line(s)
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = line[i][0].x, e = line[i][line_end[i]-1].x;
            m = i;
        } else if (L_MF(line, line_end, i) & L_MERGH) {
            r = s, l = L_RB(line, line_end, i);
            for (j = i-1; j >= m; --j) {
                if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
                E_LB(line, line_end, j) = l;
            }
            s = L_LB(line, line_end, i), e = L_RB(line, line_end, i);
            m = i;
        } //merged and body
        E_RB(line, line_end, i) = r;
    }
    for (j = i-1; j >= m; --j) {
        if (L_MF(line, line_end, j) & (L_DUMP|L_INTER)) continue;
        E_LB(line, line_end, j) = left;
    }
    // set inter-bound
    //line_inter_bound(line, line_end, li, len, t_node, tri_n, f_node, a_msg);
    return 0;
}


//for best coverage and connect
//add "overlap ins"
int get_fseed_dis(aln_msg *a_msg, int pre, int pre_a, int i, int j, int *flag, lsat_aln_per_para APP, lsat_aln_para AP)    //(i,j)对应节点，来自pre的第pre_a个aln
{
    if (pre == -1 || i == -1) { *flag = F_MATCH; return 0; }	//for bound node
    if (pre == i) {// XXX overlap???
        if (pre_a == j) *flag = F_MATCH; 
        else *flag = F_UNCONNECT;
        return 0;
    }
    if (a_msg[i].at[j].chr != a_msg[pre].at[pre_a].chr || a_msg[i].at[j].nsrand != a_msg[pre].at[pre_a].nsrand)	//different chr or different srnad
    {
        *flag = F_CHR_DIF; return 0;//PRICE_DIF_CHR;
    }

    int seed_len = APP.seed_len; int seed_inv = APP.seed_inv;
    int64_t exp = a_msg[pre].at[pre_a].offset + a_msg[pre].at[pre_a].nsrand * (a_msg[i].read_id - a_msg[pre].read_id) * (seed_len+seed_inv);	
    int64_t act = a_msg[i].at[j].offset;
    int dis = a_msg[pre].at[pre_a].nsrand * ((a_msg[pre].read_id < a_msg[i].read_id)?(act-exp):(exp-act)) - (((a_msg[pre].at[pre_a].nsrand) * (a_msg[pre].read_id-a_msg[i].read_id) < 0)?(a_msg[pre].at[pre_a].len_dif):(a_msg[i].at[j].len_dif));

    int mat_dis = APP.match_dis;
    if (dis <= mat_dis && dis >= -mat_dis) {
        if (abs(a_msg[pre].read_id - a_msg[i].read_id) == 1) *flag = F_MATCH;
        else if (abs(a_msg[pre].read_id - a_msg[i].read_id) < 10) *flag = F_MISMATCH; // XXX long dis
        else *flag = F_LONG_MISMATCH;
    } else if (dis > mat_dis && dis < AP.SV_len_thd) *flag = F_DELETE;
    else if ((dis < -mat_dis && dis >= (0-(abs(a_msg[i].read_id-a_msg[pre].read_id)*(seed_len+seed_inv)-seed_len))) // nonoverlaped ins
          || (dis < -AP.split_len && dis >= -AP.SV_len_thd)) { // overlaped ins XXX
        *flag = F_INSERT; 
    } else { *flag = F_UNCONNECT; return 0; }
    dis=abs(dis); dis += (((a_msg[pre].at[pre_a].nsrand) * (a_msg[pre].read_id - a_msg[i].read_id) < 0)? a_msg[i].at[j].cigar_len : a_msg[pre].at[pre_a].cigar_len);
    return dis; 
}

void fnode_set(frag_dp_node *f_node, line_node from, 
               int score, int dis_pen, 
               uint8_t match_flag, int dp_flag)
{
    (*f_node).son_flag = F_INIT;
    (*f_node).from = from;
    (*f_node).score = score;
    (*f_node).dis_pen = dis_pen;

    (*f_node).match_flag = match_flag;
    (*f_node).dp_flag = dp_flag;

    // general init
    (*f_node).node_n = 1;
    (*f_node).in_de = 0;
    (*f_node).son_n = 0;
    (*f_node).max_score = score;
    (*f_node).max_node.x = (*f_node).seed_i, (*f_node).max_node.y = (*f_node).aln_i;
    (*f_node).trg_n = 0;
}

int frag_dp_init(frag_dp_node **f_node, 
                 aln_msg *a_msg, int seed_i, 
                 line_node from, 
                 lsat_aln_per_para APP, lsat_aln_para AP,
                 int dp_flag)
{
    int i;
    if (from.x == START_NODE.x) { //UNLIMITED
        for (i = 0; i < a_msg[seed_i].n_aln; ++i)
            //init-score: 0 == 1 + 1 - F_SV_PEN
            fnode_set(&(f_node[seed_i][i]), from, 1, a_msg[seed_i].at[i].cigar_len, F_MATCH, dp_flag);
    } else {
        int con_flag, con_score, dis;
        for (i = 0; i < a_msg[seed_i].n_aln; ++i) {
            
            dis = get_fseed_dis(a_msg, from.x, from.y, seed_i, i, &con_flag, APP, AP);
            //											XXX MATCH, MISMATCH: same score?
            //											mismatch too long, socre??? XXX
            con_score = APP.frag_score_table[con_flag];

            if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
                fnode_set(&(f_node[seed_i][i]), from, 2+con_score, dis, con_flag, dp_flag);
            else f_node[seed_i][i].dp_flag = 0-dp_flag;
        }
    }
    return 0;
}

int fnode_add_son(frag_dp_node **f_node,
                  line_node fa, line_node son)
{
    ++f_node[fa.x][fa.y].in_de;
    int son_n = f_node[fa.x][fa.y].son_n;
    if (son_n == f_node[fa.x][fa.y].son_max) {
        f_node[fa.x][fa.y].son_max <<= 1;
        f_node[fa.x][fa.y].son = (line_node*)realloc(f_node[fa.x][fa.y].son, f_node[fa.x][fa.y].son_max * sizeof(line_node));
    }
    f_node[fa.x][fa.y].son[son_n] = son;
    ++f_node[fa.x][fa.y].son_n;
    return 0;
}

//pruning	XXX
//for best coverage and connect
//    similar to MEM, match/mismatch have the highest priority in dp.
int frag_dp_update(frag_dp_node **f_node, aln_msg *a_msg, 
                   int seed_i, int aln_i, int start, 
                   lsat_aln_per_para APP, lsat_aln_para AP, 
                   int dp_flag)
{
    int i, j, con_flag, con_score, dis;
    line_node max_from;
    int max_score, max_flag, max_dis;

    max_from = f_node[seed_i][aln_i].from;
    max_score = f_node[seed_i][aln_i].score;
    max_dis = 0; max_flag = 0;

    for (i = seed_i - 1; i >= start; --i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == dp_flag)// || f_node[i][j].dp_flag == UPDATE_FLAG)
            {
                dis = get_fseed_dis(a_msg, i, j, seed_i, aln_i, &con_flag, APP, AP);
                if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF)
                    continue;
                // '+' strand: successor-Match
                //if (con_flag > f_node[i][j].son_flag+1) //
                 //   continue;
                if (a_msg[i].at[j].nsrand == 1) {
                    if (f_node[i][j].son_flag <= F_MISMATCH) {
                        continue;
                        /// cigar-diff
                        /*if (con_flag > F_MISMATCH)
                            continue;
                        else { 
                            max_from = (line_node){i,j};
                            con_score = APP.frag_score_table[con_flag];
                            max_score = f_node[i][j].score + 1 + con_score;
                            max_flag = con_flag;
                            max_dis = f_node[i][j].dis_pen + dis;
                            goto UPDATE;
                        }*/
                    } 
                }
                // '-' strand: precursor-Match 
                ///XXX LONG_MISMATCH
                ///XXX cigar-diff, when score equal with each other
                else if (con_flag <= F_MISMATCH) {
                    max_from = (line_node){i,j};
                    con_score = APP.frag_score_table[con_flag];
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                    goto UPDATE;
                }

                con_score = APP.frag_score_table[con_flag];
                if (f_node[i][j].score + 1 + con_score > max_score)	{
                    max_from = (line_node){i,j};
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_dis = f_node[i][j].dis_pen + dis;
                } else if (f_node[i][j].score + 1 + con_score == max_score) {
                    if (f_node[i][j].dis_pen + dis < max_dis)
                    {
                        max_from = (line_node){i,j};
                        max_score = f_node[i][j].score + 1 + con_score;
                        max_flag = con_flag;
                        max_dis = f_node[i][j].dis_pen + dis;
                    }
                }
            }
        }
    }
UPDATE:
    if (max_from.x != f_node[seed_i][aln_i].from.x || max_from.y != f_node[seed_i][aln_i].from.y) {
        f_node[max_from.x][max_from.y].son_flag = max_flag;
        f_node[seed_i][aln_i].from = max_from;
        f_node[seed_i][aln_i].score = max_score;
        f_node[seed_i][aln_i].dis_pen = max_dis;
        f_node[seed_i][aln_i].match_flag = max_flag;
        f_node[seed_i][aln_i].node_n = f_node[max_from.x][max_from.y].node_n+1;

        fnode_add_son(f_node, max_from, (line_node){seed_i, aln_i}); //generate tree
    }
    return 0;
}

int frag_dp_per_init(frag_dp_node **f_node, aln_msg *a_msg, 
                     int seed_i, int aln_i, line_node from, 
                     lsat_aln_per_para APP, lsat_aln_para AP,
                     int dp_flag)
{
    if (from.x == START_NODE.x) { //UNLIMITED
        //init-score: 0 == 1 + 1 - F_SV_PEN
        fnode_set(&(f_node[seed_i][aln_i]), from, 1, a_msg[seed_i].at[aln_i].cigar_len, F_MATCH, dp_flag);
    } else {
        int con_flag, con_score, dis;
        dis = get_fseed_dis(a_msg, from.x, from.y, seed_i, aln_i, &con_flag, APP, AP);
        //											XXX MATCH, MISMATCH: same score?
        //											mismatch too long, socre??? XXX
        con_score = APP.frag_score_table[con_flag];

        if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
            fnode_set(&(f_node[seed_i][aln_i]), from, 2+con_score, dis, con_flag, dp_flag);
        else f_node[seed_i][aln_i].dp_flag = 0-dp_flag;
    }
    return 0;
}

//TODO: save all (node,score)? OR just save the first `M` big (node,score)
void node_add_score(int score, line_node node, node_score *ns)
{
    // XXX
    if (score < ns->min_score_thd) return;
    if (ns->node_n <= ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->node[(ns->node_n)++] = node;
    } else {
        fprintf(stderr, "[lsat_aln] node_add_score ERROR. (%d %d)\n", ns->node_n, ns->max_n); exit(1);
    }
}

//1.F_MATCH/F_MISMATCH/F_LONG_MISMATCH
//2.F_INS/F_DEL:
//  a.Better-Score
//  b.Short-Distance
line_node get_max_son(frag_dp_node **f_node, int x, int y)
{
    line_node max, son;
    int i, max_score = 0, max_flag, max_dis=f_node[x][y].son[0].x - x;
    for (i = 0; i < f_node[x][y].son_n; ++i) {
        son = f_node[x][y].son[i];
        max_flag = f_node[son.x][son.y].match_flag;
        /// XXX LONG_MIS_MATCH
        if (max_flag == F_MATCH || max_flag == F_MISMATCH)// || max_flag == F_LONG_MISMATCH)
            return son;
        if (f_node[son.x][son.y].max_score > max_score || (f_node[son.x][son.y].max_score == max_score && son.x - x < max_dis)) {
            max = son;
            max_score = f_node[son.x][son.y].max_score;
            max_dis = son.x - x;
        } 
        /*if (son.x - x == max_dis) { // same dis, better score
            if (f_node[son.x][son.y].max_score > max_score)
            {
                max = son;
                max_score = f_node[son.x][son.y].max_score;
            }
        } else if (son.x - x < max_dis) {
            max = son;
            max_score = f_node[son.x][son.y].max_score;
            max_dis = son.x - x;
        }*/
    }
    return max;
}

void cut_branch(frag_dp_node **f_node, int x, int y, node_score *ns) {
    line_node max_son = get_max_son(f_node, x, y);
    line_node son, max_node;
    int i;
    for (i = 0; i < f_node[x][y].son_n; ++i) {
        son = f_node[x][y].son[i];
        if (nodeEq(max_son, son)) continue;
        //update score: 
        //New_Max = Old_Max - branch_score - con_score 
        //        = Old_Max - son_score + 1
        //Because: son_score = branch_score + 1 + con_score
        f_node[son.x][son.y].from = START_NODE;
        f_node[son.x][son.y].max_score -= (f_node[son.x][son.y].score - 1);
        max_node = f_node[son.x][son.y].max_node;
        f_node[max_node.x][max_node.y].node_n -= (f_node[son.x][son.y].node_n-1);
        node_add_score(f_node[son.x][son.y].max_score, max_node , ns);
    }
    // refresh max_socre/node of (x,y), based on max_son, consider negative edge
    if (f_node[x][y].score > f_node[max_son.x][max_son.y].max_score) { // negative edge
        f_node[max_son.x][max_son.y].in_de = -1;
        f_node[max_son.x][max_son.y].from = START_NODE;
        f_node[max_son.x][max_son.y].max_score -= (f_node[max_son.x][max_son.y].score - 1);
        max_node = f_node[max_son.x][max_son.y].max_node;
        f_node[max_node.x][max_node.y].node_n -= (f_node[max_son.x][max_son.y].node_n-1);
        node_add_score(f_node[max_son.x][max_son.y].max_score, max_node, ns);
        f_node[x][y].son_n = 0;
        f_node[x][y].max_node = (line_node){x,y};
        f_node[x][y].max_score = f_node[x][y].score;
    } else {
        f_node[x][y].son_n = 1;
        f_node[x][y].son[0] = max_son;
        f_node[x][y].max_node = f_node[max_son.x][max_son.y].max_node;
        f_node[x][y].max_score = f_node[max_son.x][max_son.y].max_score;
    }
    f_node[x][y].in_de = 0; // leaf-node, avaible for branch-track
}

//trackback from (x,y)
void branch_track_new(frag_dp_node **f_node, int x, int y, node_score *ns)
{
    int max_score;line_node max_node; 
    f_node[x][y].in_de = -1;
    if (f_node[x][y].son_n == 0) {
        max_node = f_node[x][y].max_node = (line_node){x,y};
        max_score = f_node[x][y].max_score = f_node[x][y].score;
    } else {
        max_node = f_node[x][y].max_node;
        max_score = f_node[x][y].max_score;
    }

    line_node fa = f_node[x][y].from, son; int fa_son_n;

    while (fa.x != START_NODE.x) {
        fa_son_n = f_node[fa.x][fa.y].son_n;
        if (fa_son_n == 1) { // 继续回溯
            if (f_node[fa.x][fa.y].score > max_score) { //negative edge
                son = f_node[fa.x][fa.y].son[0];
                f_node[son.x][son.y].in_de = -1;
                f_node[son.x][son.y].from = START_NODE;
                f_node[son.x][son.y].max_score -= (f_node[son.x][son.y].score - 1);
                f_node[max_node.x][max_node.y].node_n -= (f_node[son.x][son.y].node_n-1);
                node_add_score(f_node[son.x][son.y].max_score, max_node, ns);
                f_node[fa.x][fa.y].son_n = 0;
                max_score = f_node[fa.x][fa.y].score;
                max_node = fa;
            }
            f_node[fa.x][fa.y].max_score = max_score;
            f_node[fa.x][fa.y].max_node = max_node;
            f_node[fa.x][fa.y].in_de = -1;
            fa = f_node[fa.x][fa.y].from;
        } else { // fa_son_n > 1, TREE_BRANCH
            --f_node[fa.x][fa.y].in_de;
            if (f_node[fa.x][fa.y].in_de == 0)
                cut_branch(f_node, fa.x, fa.y, ns);
            return;
        }
    }
    // TREE_START; 
    node_add_score(max_score, max_node, ns);// add (max_node, max_score);
    return;
}



//for multi-dp-lines
//line, line_end: start at 1
//left and right are 'MIN_FLAG'	XXX
int frag_mini_dp_multi_line(frag_dp_node **f_node, aln_msg *a_msg, 
                            lsat_aln_per_para APP, lsat_aln_para AP,
                            int left_b, int right_b, 
                            line_node **line, int *line_end,
							int line_n_max)
{
    if (left_b+1 >= right_b) return 0; //XXX 

    line_node head = START_NODE; // unlimited
    int start, end;
    int i, j, mini_dp_flag = WHOLE_FLAG;// (whole DP)
    int l_i, node_i;
    line_node left, right;
    line_node _right;
    //int candi_n; line_node candi[max_multi];

    left.x = left_b, right.x = right_b;
    start = left.x+1;
    end = right.x-1;

    //first dp int
    for (i = start; i <= end; ++i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            //if (f_node[i][j].dp_flag == mini_dp_flag || f_node[i][j].dp_flag == 0-mini_dp_flag)
                frag_dp_per_init(f_node, a_msg, i, j, head, APP, AP, mini_dp_flag);
        }
    }
    //dp update
    for (i = start+1; i <= end; ++i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag)
                frag_dp_update(f_node, a_msg, i, j, start, APP, AP, mini_dp_flag);
        }
    }
    // tree-pruning
    node_score *ns = node_init_score(line_n_max); //size of node_score XXX
    ns->min_score_thd = 0;
    //node_score *ns = node_init_score(80); //size of node_score
    for (i = end; i >=start; --i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag && f_node[i][j].in_de == 0) //leaf node
                branch_track_new(f_node, i, j, ns); 
        }
    }
    l_i = 0;
    int score;
    while (1) {
        _right = node_pop(ns, &score);
        if (_right.x == START_NODE.x) break;
        node_i = f_node[_right.x][_right.y].node_n-1;
        line_end[l_i] = node_i+1;
        L_LS(line,line_end,l_i) = score;
        L_BS(line,line_end,l_i) = score;
        while (_right.x != head.x) {
            if (node_i < 0) { 
                fprintf(stderr, "\n[frag mini dp multi] %s node_i BUG 1.\n", APP.read_name); exit(1); 
			}
            line[l_i][node_i--] = _right;
            f_node[_right.x][_right.y].dp_flag = MULTI_OUT;
            _right = f_node[_right.x][_right.y].from;
        }
        if (node_i >= 0) { 
            fprintf(stderr, "\n[frag mini dp multi] %s node_i BUG 2.\n", APP.read_name); exit(1); 
		}
        ++l_i;
    }
    node_free_score(ns);
    return l_i;
}

int trg_dp_line(frag_dp_node **f_node, aln_msg *a_msg,
                lsat_aln_per_para APP, lsat_aln_para AP, 
                int left, int right, 
                line_node **line, int *line_end, 
                int line_n_max, int per_max_multi) 
{
    int l = frag_mini_dp_multi_line(f_node, a_msg, APP, AP, left, right, line, line_end, line_n_max);//, 0, 0);
    line_set_bound1(line, line_end, 0,  &l, left, right, per_max_multi);
    ///l = line_remove(line, line_end, 0, l);
    return l;
}

void frag_min_extend(frag_dp_node **f_node, aln_msg *a_msg,
                    int node_i, int aln_i,
                    int aln_min, int dp_flag,
                    lsat_aln_per_para APP, lsat_aln_para AP)
{
    int i, j, con_flag;
    int last_x = node_i, last_y = aln_i;
    // 0 -> node_i-1, node_i+1 -> node_n-1
    //from right to left 
    i = last_x - 1;
    while (i >= 0) {
        if (a_msg[i].n_aln > aln_min) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                get_fseed_dis(a_msg, i, j, last_x, last_y, &con_flag, APP, AP);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH || con_flag == F_LONG_MISMATCH) {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        --i;
    }
    i = last_x + 1;
    while (i < APP.seed_out) {
        if (a_msg[i].n_aln > aln_min) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                get_fseed_dis(a_msg, last_x, last_y, i, j, &con_flag, APP, AP);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH || con_flag == F_LONG_MISMATCH) {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        ++i;
    }
}

//mini-dp 过程中，种子比对结果被多重输出的问题
//解决办法:每次mini-dp结束后，修改该次输出的node对应的dp-flag，
//         以保证在后面的mini-dp中，它们不再被使用到.
//for best coverage and connect
//fill-up for inter-gaps in min-line: XXX 如何填补? 如何避免出现错误填补,即将不该填补的补上了
int frag_mini_dp_line(frag_dp_node **f_node, aln_msg *a_msg, 
                      lsat_aln_per_para APP, lsat_aln_para AP,
                      line_node left, line_node right, 
                      line_node *line, int *de_score,
                      int _head, int _tail)
{
    extern int f_BCC_score_table[10];
    line_node head = ((_head)?left:START_NODE);
    int old_score;
    if (_head==0 || _tail==0) old_score = 1;
    else old_score = 2+f_BCC_score_table[f_node[right.x][right.y].match_flag];
    int i, j, mini_dp_flag = MULTI_FLAG;
    //dp init
    for (i = left.x + 1; i < right.x; ++i) {
        //mini_dp XXX?
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag || f_node[i][j].dp_flag == 0-mini_dp_flag)
                frag_dp_per_init(f_node, a_msg, i, j, head, APP, AP, mini_dp_flag);
        }
    }
    //dp update
    for (i = left.x + 2; i < right.x; ++i) {
        for (j = 0; j < a_msg[i].n_aln; ++j) {
            if (f_node[i][j].dp_flag == mini_dp_flag)
                frag_dp_update(f_node, a_msg, i, j, left.x+1, APP, AP, mini_dp_flag);
        }
    }
    //find backtrack start node
    int max_score, max_dis = 0, max_n = 0;
    line_node max_node = head;//left;
    //UNLIMITED 
    if (_tail == 0) {
        max_score = old_score;
        for (i = right.x - 1; i > left.x; --i) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                if (f_node[i][j].dp_flag == mini_dp_flag) {
                    if (f_node[i][j].score > max_score) {
                        max_score = f_node[i][j].score;
                        max_dis = f_node[i][j].dis_pen;
                        max_node = (line_node){i, j};
                        max_n = f_node[i][j].node_n;
                    } else if (f_node[i][j].score == max_score) {
                        if (f_node[i][j].dis_pen < max_dis) {
                            max_score = f_node[i][j].score;
                            max_dis = f_node[i][j].dis_pen;
                            max_node = (line_node){i, j};
                            max_n = f_node[i][j].node_n;
                        }
                    }
                }
            }
        }
    } else {
        f_node[right.x][right.y].from = head;//left;
        f_node[right.x][right.y].score = old_score;
        f_node[right.x][right.y].node_n = 1;
        frag_dp_update(f_node, a_msg, right.x, right.y, left.x+1, APP, AP, mini_dp_flag);
        max_score = f_node[right.x][right.y].score;
        max_node = f_node[right.x][right.y].from;
        max_n = f_node[right.x][right.y].node_n-1; // right is NOT inculding 
    }
    //backtrack
    line_node _right = max_node;
    int node_i = max_n-1;
    //while (_right.x != left.x)
    while (_right.x != head.x) {
        if (node_i < 0) { 
            fprintf(stderr, "\n[frag mini dp] %s node_i BUG 1.\n", APP.read_name); exit(1); 
        }
        line[node_i--] = _right;
        _right = f_node[_right.x][_right.y].from;
    }
    if (node_i >= 0) { fprintf(stderr, "\n[frag mini dp] %s node_i BUG 2.\n", APP.read_name); exit(1); }
    *de_score += (max_score-old_score);
    return max_n;
}

int frag_dp_path(aln_msg *a_msg, frag_msg **f_msg,
                 lsat_aln_per_para APP,
                 int *l_n, int *line_m,
                 line_node **line, int *line_end,
                 frag_dp_node ***f_node)
{
    int line_n = *l_n;
    if (line_n == 0) return 0;
	int i, l;
    
#ifdef __DEBUG__
    for (i = 0; i < line_n; ++i) {
        printf("%d(%d,%d score: %d):\t", i+1, L_MF(line, line_end, i), L_MH(line, line_end,i), L_LS(line, line_end, i));
        for (l = 0; l < line_end[i]; ++l)
            printf("(%d, %d)\t", a_msg[line[i][l].x].read_id-1, line[i][l].y);
            //printf("(%d, %d)\t", line[i][l].x, line[i][l].y);
        printf("\n");
    }
#endif

	int frag_num;
	int cur_x , cur_y , pre_x, pre_y;

	int left_bound, right_bound;

    if (line_n > (*line_m)) {
        (*f_msg) = (frag_msg*)realloc(*f_msg, line_n * sizeof(frag_msg));
        for (i = (*line_m); i < line_n; ++i)
        {
            frag_init_msg((*f_msg)+i, (*f_msg)->frag_max); 
            frag_copy_msg(*f_msg, (*f_msg)+i);
        }
        if ((*f_msg) == NULL) { fprintf(stderr, "\n[frag_dp_path] Not enough memory.(line_m: %d)\n", line_n); exit(1); }
        (*line_m)= line_n;
        fprintf(stderr, "line-num: %d\t", line_n);
    }

    int j=0;
	for (l = 0; l < line_n; ++l) {
        //set merged msg
        ((*f_msg)+j)->merg_msg.x = L_MF(line,line_end,l);
        ((*f_msg)+j)->merg_msg.y = L_MH(line,line_end,l);
            
        frag_num=0;
        pre_x = line[l][line_end[l]-1].x; pre_y = line[l][line_end[l]-1].y;
        //fprintf(stdout, "left: %d, right: %d\n", line[l][line_end[l]].x, line[l][line_end[l]+1].x);
		//right bound
        //right_bound = ((E_RB(line, line_end, l) == APP.seed_out) ? (APP.seed_all+1) : a_msg[E_RB(line, line_end, l)].read_id);
		right_bound = APP.seed_all+1;
        //MIS-MATCH
        //first end
		frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num); //XXX set merged msg
		if ((*f_node)[pre_x][pre_y].trg_n)	//set trigger
            frag_trg_set((*f_node)[pre_x][pre_y], (*f_msg)+j, frag_num);
        ((*f_msg)+j)->frag_right_bound = right_bound;
		for (i = line_end[l]-1; i > 0; --i) {
			cur_x = pre_x; cur_y = pre_y;
			pre_x = line[l][i-1].x; pre_y = line[l][i-1].y;
			//INS
			if ((*f_node)[cur_x][cur_y].match_flag == F_INSERT)
			{
                frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
				if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
                    frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num);
            }
			//DEL
			else if ((*f_node)[cur_x][cur_y].match_flag == F_DELETE) {
				frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
				if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
                    frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num);
			}
			//MIS	XXX
			else if ((*f_node)[cur_x][cur_y].match_flag == F_MISMATCH || (*f_node)[cur_x][cur_y].match_flag == F_LONG_MISMATCH) {
				//XXX take mis-match as NEW-FLAG case
				//frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
				frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
				if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
                    frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
                ++frag_num;
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_END, (*f_msg)+j, frag_num);
				//XXX find new flag in mis-match seeds
				//dp-line for INV/TRS
            } else if ((*f_node)[cur_x][cur_y].match_flag == F_MATCH)	//: MATCH
			{
				frag_set_msg(a_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+j, frag_num);
			} else { fprintf(stderr, "\n[frag dp path] Error: Unknown flag, \"%d\"\n", (*f_node)[cur_x][cur_y].match_flag); exit(1);}
		}
		//last start
		cur_x = line[l][0].x; cur_y = line[l][0].y;
		frag_set_msg(a_msg, cur_x, cur_y, FRAG_START, (*f_msg)+j, frag_num);
        if ((*f_node)[cur_x][cur_y].trg_n)	//set trigger
            frag_trg_set((*f_node)[cur_x][cur_y], (*f_msg)+j, frag_num);
        //left_bound = ((E_LB(line, line_end, l) == -1) ? 0 : a_msg[E_LB(line, line_end, l)].read_id);
		left_bound = 0;

		//MIS-MATCH
        ((*f_msg)+j)->frag_left_bound = left_bound;
#ifdef __DEBUG__
        fprintf(stdout, "#%d: left %d, right %d\n", j+1, left_bound, right_bound);
#endif
        ++j;
	}
    *l_n = j;
	
	/*int k, s_i, a_i;
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
    return 1;
}

int frag_line_remain(aln_reg *a_reg, aln_msg *a_msg, frag_msg **f_msg,
                     lsat_aln_per_para APP, lsat_aln_para AP, 
                     line_node **line, int *line_end, int *line_m, 
                     frag_dp_node ***f_node, line_node **_line, int *_line_end,
                     int line_n_max, int per_max_multi)
{
    int l_n = 0;
    aln_reg *re_reg = aln_init_reg(APP.read_len);
    if (get_remain_reg(a_reg, re_reg, AP, APP.read_len) == 0) goto End;
    int i, j, k;
    int left_id, right_id, left, right, l;
    for (i = 0; i < re_reg->reg_n; ++i) {
        // get region boundary (left, right)
        left_id = (re_reg->reg[i].beg +  APP.seed_inv - 1) / APP.seed_step+1; 
        right_id = (re_reg->reg[i].end - 1) / APP.seed_step+1;
		if (right_id > APP.seed_all) right_id -= 1;
		left = right = -2;
		for (j = 0; j < APP.seed_out; ++j) {
			if (a_msg[j].read_id >= left_id) {
				left = j-1;
				break;
			}
		}
		if (left == -2) continue;
		for (j = APP.seed_out-1; j >= 0; --j) {
			if (a_msg[j].read_id <= right_id) {
				right = j+1;
				break;
			}
		}
		if (right == -2) continue;
		// store line-info in _line
		l = trg_dp_line(*f_node, a_msg, APP, AP, left, right, _line, _line_end, line_n_max, per_max_multi);
		// copy from _line to line
        for (j = 0; j < l; ++j) {
            line_end[l_n+j] = _line_end[j];
            for (k = 0; k < _line_end[j]+L_EXTRA; ++k) {
                line[l_n+j][k] = _line[j][k];
            }
        }
        l_n += l;
    }
End:
    aln_free_reg(re_reg);
    frag_dp_path(a_msg, f_msg, APP, &l_n, line_m, line, line_end, f_node);
    return l_n;
}

//new tree-generating and pruning
int frag_line_BCC(aln_msg *a_msg, frag_msg **f_msg,
                  lsat_aln_per_para APP, lsat_aln_para AP,
                  line_node **line, int *line_end, int *line_m,
                  frag_dp_node ***f_node,
                  line_node **_line, int line_n_max, int per_max_multi)
{
    int i, j, k, line_n;
    int min_n = APP.min_thd, min_exist=0, min_num = 0;

    { // dp init
        for (i = 0; i < APP.seed_out; ++i)
        {
            if (a_msg[i].n_aln <= min_n)
            {
                frag_dp_init(*f_node, a_msg, i, START_NODE, APP, AP, MIN_FLAG);
                min_exist = 1;
                ++min_num;
            } else frag_dp_init(*f_node, a_msg, i, START_NODE, APP, AP, MULTI_FLAG);
        }
        //fraction XXX
        if (!min_exist || min_num * 3 < APP.seed_out)
        {
            for (i = 0; i < APP.seed_out; ++i) {
                for (j = 0; j < a_msg[i].n_aln; ++j) (*f_node)[i][j].dp_flag = MIN_FLAG;
            }
            min_n = APP.per_aln_n;
            min_exist = 1;
        }
    }
    { // dp update
        //min extend, when min_n == PER_ALN_N: no need to extend
        if (min_n != APP.per_aln_n)
        {
            for (i = 0; i < APP.seed_out; ++i) {
                if (a_msg[i].n_aln <= min_n) {
                    for (j = 0; j < a_msg[i].n_aln; ++j)
                        frag_min_extend(*f_node, a_msg, i, j, min_n, MIN_FLAG, APP, AP);
                }
            }
        }
        //min update
        for (i = 1; i < APP.seed_out; ++i) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                if ((*f_node)[i][j].dp_flag == MIN_FLAG)
                    frag_dp_update(*f_node, a_msg, i, j, 0/*update start pos*/, APP, AP, MIN_FLAG);
            }
        }
    }
    { // backtrack
        //prune tree to find multi-nodes/roads for backtrack
        //XXX
        node_score *ns = node_init_score(line_n_max); //size of node_score
        ns->min_score_thd = 2;
        for (i = APP.seed_out-1; i >=0; --i) {
            for (j = 0; j < a_msg[i].n_aln; ++j) {
                if ((*f_node)[i][j].dp_flag == MIN_FLAG && (*f_node)[i][j].in_de == 0) //leaf node
                    branch_track_new(*f_node, i, j, ns); 
            }
        }
        int line_score; line_node max_node;
        int l_i = 0, new_l = 0, min_l=ns->node_n, o_l = min_l;
        int node_i, mini_len;
        trig_node **inter_tri = (trig_node**)malloc(o_l* sizeof(trig_node*)); int *inter_n=(int*)malloc(o_l*sizeof(int));
        for (i = 0; i < o_l;++i)
            inter_tri[i] = (trig_node*)malloc(APP.seed_out * sizeof(trig_node));
        //multi-backtrack
        //extend for min-line AND set inter-trigger
        while (1) {
            //MAX score //XXX
            max_node = node_pop(ns, &line_score);    //这里得到的path,就是希望保留的
                                                     //总共ns->node_n条path,即min-line
            if (max_node.x == -1) break;
            //backtrack for max-node
            line_node last_n, right, left;
            //fill-up for right end
            node_i=0;
            inter_n[l_i]=0;
            if (max_node.x < APP.seed_out-1) {
                //XXX 是否有意义?
                mini_len = frag_mini_dp_line(*f_node, a_msg, APP, AP, max_node, (line_node){APP.seed_out, 0}, _line[0], &line_score, 1, 0);
                for (k = mini_len-1; k >= 0; --k) {
                    line[l_i][node_i++] = _line[0][k];
                    (*f_node)[_line[0][k].x][_line[0][k].y].dp_flag = MULTI_OUT; //不允许重复使用
                }
                line[l_i][node_i] = max_node;
                //set inter-trigger
                last_n = line[l_i][0];
                for (k = mini_len-1; k >= 0; --k) {
                    if (last_n.x - line[l_i][node_i-k].x > 2) {
                        // pre-trigger: (line[l_i][node_i-k].x , right.x)
                        (*f_node)[last_n.x][last_n.y].pre_trg = (line_node){line[l_i][node_i-k].x, last_n.x};
                        // 同一node处有多个trigger?XXX,即trg_n被重复赋值1或2
                        (*f_node)[last_n.x][last_n.y].trg_n |= 0x2;
                        inter_tri[l_i][inter_n[l_i]++] = (trig_node){(line_node){line[l_i][node_i-k].x, line[l_i][node_i-k].y}, (line_node){last_n.x, last_n.y}, 2};
                    }
                    last_n = line[l_i][node_i-k];
                }
            }
            right = max_node;
            // fill-up for every gap
            while (right.x != START_NODE.x) {
                line[l_i][node_i++] = right;
                left = (*f_node)[right.x][right.y].from;
                if (left.x < right.x - 1 )//f_node[right.x][right.y].match_flag != F_MATCH)	//XXX if left.x < right.x-1, match_flag of right node couldn't be match?
                {
                    //mini-dp with anchors
                    mini_len = frag_mini_dp_line(*f_node, a_msg, APP, AP, left, right, _line[0], &line_score, 1, 1);
                    for (k = mini_len-1; k >= 0; --k) {
                        line[l_i][node_i++] = _line[0][k];
                        (*f_node)[_line[0][k].x][_line[0][k].y].dp_flag = MULTI_OUT; //不允许重复使用
                    }
                    line[l_i][node_i] = left;	//twice write, for the last multi-dp-line	//XXX
                    //set inter-trigger
                    last_n = right;
                    for (k = mini_len; k >= 0; --k) {
                        if (last_n.x - line[l_i][node_i-k].x > 2) {
                            // pre-trigger: (line[l_i][node_i-k].x , right.x)
                            if (line[l_i][node_i-k].x == START_NODE.x) continue;
                            (*f_node)[last_n.x][last_n.y].pre_trg = (line_node){line[l_i][node_i-k].x, last_n.x};
                            (*f_node)[last_n.x][last_n.y].trg_n |= 0x2;
                            inter_tri[l_i][inter_n[l_i]++] = (trig_node){(line_node){line[l_i][node_i-k].x, line[l_i][node_i-k].y},(line_node){last_n.x, last_n.y}, 2};
                        }
                        last_n = line[l_i][node_i-k];
                    }
                }
                right = left;
            }
            //invert line
            line_node tmp;
            for (k = 0; k < node_i/2; ++k) {
                tmp =line[l_i][k]; line[l_i][k] = line[l_i][node_i-k-1]; line[l_i][node_i-k-1] = tmp;
            }
            line_end[l_i] = node_i;
            //line[l_i][line_end[l_i]+2].x = line_score;
            L_LS(line, line_end, l_i) = line_score; L_BS(line, line_end, l_i) = line_score;
            l_i++;
        }
        /* line: [0 ~ end-1] -> node; 
                 [end] -> (left-b, right-b)/(start,end)
                 [end+1] -> (merged-flag, merge-head)/(start, end) // (-1,):NOT merged/ (-2,i):Merged and head is i/ (a,b): start from a to b
                 [end+2] -> (line-score, x) */
        // min-lines have been extended, and sorted by end-pos
        // set boundary and filter best/secondary
        line_set_bound(line, line_end, 0, &min_l, -1, APP.seed_out, inter_tri, inter_n, *f_node, a_msg, per_max_multi);
        for (i = 0; i < o_l; ++i) free(inter_tri[i]);
        free(inter_tri); free(inter_n);
        ///min_l = line_remove(line, line_end, 0, min_l);

        // set bound-trigger (-1, right)/ (left, APP.seed_out)
        for (i = 0; i < min_l; ++i) {
            if(L_MF(line, line_end, i) & L_INTER) continue;
            if (E_RB(line, line_end, i) == APP.seed_out && line[i][line_end[i]-1].x < APP.seed_out-2) {
                (*f_node)[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].next_trg = (line_node){line[i][line_end[i]-1].x, APP.seed_out};
                (*f_node)[line[i][line_end[i]-1].x][line[i][line_end[i]-1].y].trg_n |= 0x1;
            } 
            if (E_LB(line, line_end, i) == -1 && line[i][0].x > 1) {
                (*f_node)[line[i][0].x][line[i][0].y].pre_trg = (line_node){-1, line[i][0].x};
                (*f_node)[line[i][0].x][line[i][0].y].trg_n |= 0x2;
            }
        }
        node_free_score(ns);
        line_n = min_l+new_l;

        frag_dp_path(a_msg, f_msg, APP, &line_n, line_m, line, line_end, f_node);
        return line_n;
    }
}
