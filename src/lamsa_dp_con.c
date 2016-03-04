#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lamsa_aln.h"
#include "lamsa_dp_con.h"
#include "split_mapping.h"
#include "lamsa_heap.h"

typedef struct { int pos; int li; } line_pos_t;
int line_pos_comp(const void *a, const void *b) { return ((line_pos_t*)b)->pos - ((line_pos_t*)a)->pos; }
void line_sort_endpos(line_node *line, int *lsl, 
                      int *line_rank, int *line_select_rank,
                      int ls, int len) {
    int i;
    line_pos_t *pt = (line_pos_t*)malloc(len * sizeof(line_pos_t));

    for (i = ls; i < ls+len; ++i) {
        pt[i-ls].li = i, pt[i-ls].pos = _line_node(line, lsl, i)[_line_len(lsl, i)-1].x;
    }
    qsort(pt, len, sizeof(line_pos_t), line_pos_comp);
    for (i = ls; i < ls+len; ++i) {
        line_rank[i] = pt[i-ls].li;
        line_select_rank[pt[i-ls].li] = i;
    }
    free(pt);
}

node_score *node_init_score(int n) {
    node_score *ns = (node_score*)malloc(sizeof(node_score));
    ns->max_n = n;
    ns->node_n = 0;
    ns->node = (line_node*)malloc(n * sizeof(line_node));
    ns->score = (int*)malloc(n * sizeof(int));
    ns->NM = (int*)malloc(n * sizeof(int));
    return ns;
}

void node_free_score(node_score *ns) {
    free(ns->score); free(ns->NM); free(ns->node);
    free(ns);
}

int heap_add_node(node_score *ns, line_node node, int score, int NM) {
    if (ns->node_n < ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->NM[ns->node_n] = NM;
        ns->node[ns->node_n++] = node;
        return -1;
    } else if (ns->node_n == ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->NM[ns->node_n] = NM;
        ns->node[ns->node_n++] = node;
        build_node_min_heap(ns);
        return -1;
    } else { // (ns->node_n == max_n, min-heap is already built.
        return node_heap_update_min(ns, node, score, NM); 
    }
}

float cover_rate(int s1, int e1, int s2, int e2) {
    int s = (s2 > s1) ? s2 : s1;
    int e = (e2 < e1) ? e2 : e1;
    float rat1 = (e-s+1+0.0)/(e1-s1+1+0.0);
    float rat2 = (e-s+1+0.0)/(e2-s2+1+0.0);
    return rat1>rat2?rat1:rat2;
}

int line_merge(int a, int b, line_node *line, int *lsl, float ovlp_r) {
    int s1, e1, s2, e2, s, e, hi;
    float rat1, rat2;
    line_node *na = _line_node(line, lsl, a), *nb = _line_node(line, lsl, b);
    int la = _line_len(lsl, a), lb = _line_len(lsl, b);
    line_node *nhi; int lhi;


    s2 = na[0].x, e2 = na[la-1].x;
    if (L_MF(nb,lb) & L_NMERG) {// Not merged
        hi = b;
        nhi = nb, lhi = lb;
        s1 = nb[0].x, e1 = nb[lb-1].x;
    } else { // [b] was merged already
        hi = L_MH(nb, lb);
        nhi = _line_node(line, lsl, hi), lhi = _line_len(lsl, hi);
        s1 = L_LB(nhi, lhi), e1 = L_RB(nhi, lhi);
    }
    s = (s2 > s1) ? s2 : s1;
    e = (e2 < e1) ? e2 : e1;
    
    rat1 = (e-s+1+0.0)/(e1-s1+1+0.0);
    rat2 = (e-s+1+0.0)/(e2-s2+1+0.0);
	if (rat1 < ovlp_r && rat2 < ovlp_r) {
        L_MF(na, la) = L_NMERG;
        return 0;
	} else if 
		(L_LS(na, la) <= L_LS(nb, lb) / 2 || // merged and  dumped
		 L_LS(na, la) <= L_BS(nb, lb) / 2 ) {

        L_LB(nhi, lhi) = s1, L_RB(nhi, lhi) = e1;
        L_MF(nhi, lhi) = L_MERGH;
        L_MF(na, la) = L_MERGB, L_MH(na, la) = hi;
		L_MF(na, la) |= L_DUMP;
		return 1;
	} else { // merged
        L_LB(nhi, lhi) = s1+s2-s, L_RB(nhi, lhi) = e1+e2-e;
        L_MF(nhi, lhi) = L_MERGH;
        L_MF(na, la) = L_MERGB, L_MH(na, la) = hi;
		if (L_BS(nb, lb) > L_BS(na, la))
			L_BS(na, la) =  L_BS(nb, lb);
        return 1;
    } 
}

// rule: 1. Better Score
//       2. Short line, if Score is equal
//            
// Best & Secondary: 1. Best is unique, 
//                   2. Secondary is allowed to be multi, 
//            		 3. Secondary is allowed to equal the Best
// new
// Best and all secondary whose score exceed best/2
void line_filter(line_node *line, int *lsl, int *line_rank, int *line_select_rank, int ls, int len, 
                 trig_node **trg_node, int *tri_n, frag_dp_node **f_node, 
                 map_msg *m_msg, int per_max_multi)
{
    int i, _i, ii, j, _j, k, l, _l, m_ii;
    line_node *ni; int li;
    tri_node **m_b = (tri_node**)malloc(len * sizeof(tri_node*)); // (i, score, NM)
	int **m_f = (int**)malloc(len * sizeof(int*)); // merged line (after filtering)
	for (i = 0; i < len; ++i) {
		m_b[i] = (tri_node*)malloc(len * sizeof(tri_node));
		m_f[i] = (int*)malloc((len+1) * sizeof(int));
	}
	int *m_bn = (int*)malloc(len * sizeof(int)); 
	int *m_fn = (int*)malloc(len * sizeof(int));
    int m_i = -1;
    int b_score, s_score;
	int m_head, f_merge, min_l, max_r, head;
    // check merged-lines
    for (_i = ls; _i < ls+len; ++_i) {
        i = line_rank[_i]; ni = _line_node(line, lsl, i), li = _line_len(lsl, i);
		if (L_MF(ni, li) & L_DUMP) continue;
        if (L_MF(ni, li) & L_NMERG) { // NOT merged
			++m_i;
			m_b[m_i][0].x = i, m_b[m_i][0].y = -2;
			m_bn[m_i] = 1;
		} else if (L_MF(ni, li) & L_MERGH) { // merged, head
			++m_i;
            m_b[m_i][0].x = i, m_b[m_i][0].y = L_LS(ni, li); // NO secondary
            m_b[m_i][0].z = L_NM(ni, li);
			m_bn[m_i] = 1;
		} else { // merged, body : inter-line candidate
            //if (m_i < 0) { fprintf(stderr, "[line_filter] %s BUG.\n", READ_NAME); exit(1); }
			m_b[m_i][m_bn[m_i]].x = i, m_b[m_i][m_bn[m_i]].y = L_LS(ni, li), m_b[m_i][m_bn[m_i]].z = L_NM(ni, li);
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
		if (s_score >= b_score/2) f_merge = 1; // get M largest scores
		else f_merge = 0;                     // get the largest score

		m_fn[i] = 1; // m_f [0: best, 1: head, 1...: body
		if (f_merge) {
            //1. init_heap
            ns->node_n = 0;
			for (j = 0; j < m_bn[i]; ++j) {
                //5. output_M_nodes
				if (m_b[i][j].y >= b_score/2) {
                    //2. add_node & build_heap & update_heap
                    int ret = heap_add_node(ns, (line_node){m_b[i][j].x,-1}, m_b[i][j].y, m_b[i][j].z);
                    if (ret == -2) { // dumped
                        _L_MF(line, lsl, m_b[i][j].x) |= L_DUMP;
                        tri_n[m_b[i][j].x] = 0;
                    } else if (ret != -1) { // ret is dumped
                        _L_MF(line, lsl, ret) |= L_DUMP;
                        tri_n[ret] = 0;
                    }
                } else {
					_L_MF(line, lsl, m_b[i][j].x) |= L_DUMP;
					tri_n[m_b[i][j].x] = 0;
				}
			} 	
            //3. output_M_nodes, set head/body flags
            build_node_minpos_heap(ns);
            // head
            line_node head_n = node_heap_extract_minpos(ns);
            m_head = head_n.x;
            line_node *mh_n = _line_node(line, lsl, m_head); int mh_l = _line_len(lsl, m_head);
            L_MF(mh_n, mh_l) = L_MERGH;

            min_l = mh_n[0].x;
            max_r = mh_n[mh_l-1].x;
            if (L_LS(mh_n, mh_l) == b_score) m_f[i][0] = m_head; // best
            m_f[i][m_fn[i]++] = m_head;

            int body;
            while ((body = node_heap_extract_minpos(ns).x) != -1) {
                line_node *bd_n = _line_node(line, lsl, body); int bd_l = _line_len(lsl, body);
                L_MF(bd_n, bd_l) = L_MERGB;
                L_MH(bd_n, bd_l) = m_head;
                min_l = MINOFTWO(min_l, bd_n[0].x);
                max_r = MAXOFTWO(max_r, bd_n[bd_l-1].x);
                if (L_LS(bd_n, bd_l) == b_score) m_f[i][0] = body; // best
                m_f[i][m_fn[i]++] = body;
            }
            // set left/right boundary
			L_LB(mh_n, mh_l) = MINOFTWO(mh_n[0].x, min_l);
			L_RB(mh_n, mh_l) = MAXOFTWO(mh_n[mh_l-1].x, max_r);
		} else {
			for (j = 0; j < m_bn[i]; ++j) {
				if (m_b[i][j].y == b_score) {
					_L_MF(line, lsl, m_b[i][j].x) =  L_NMERG;
					m_f[i][0] = m_b[i][j].x; // best
					m_f[i][m_fn[i]++] = m_b[i][j].x;
				} else {
					_L_MF(line, lsl, m_b[i][j].x) |= L_DUMP;
					tri_n[m_b[i][j].x] = 0;
				}
			}
		}
		for (ii = 1; ii < m_fn[i]; ++ii) {
			j = m_f[i][ii]; _j = line_select_rank[j];
			for (k = 0; k < tri_n[j]; ++k) {
                head = -1;
				// check inter-line: from j+1 -> merge_end
                for (_l = _j+1; _l < ls+len ; ++_l) { // L_MERGB
                    l = line_rank[_l];
                    line_node *nl = _line_node(line, lsl, l); int ll = _line_len(lsl, l);
                    if ((L_MF(nl, ll) & 0x3) != 0) break;
                    if (nl[0].x > trg_node[j][k].n1.x && nl[ll-1].x < trg_node[j][k].n2.x) {
                        // check for MIS(INV)
						if (f_node[trg_node[j][k].n2.x][trg_node[j][k].n2.y].match_flag == F_MISMATCH ||
                                f_node[trg_node[j][k].n2.x][trg_node[j][k].n2.y].match_flag == F_LONG_MISMATCH) {
                            int x_s = nl[0].x, y_s = nl[0].y;
                            int x_e = nl[ll-1].x, y_e = nl[ll-1].y;
                            int x1 = trg_node[j][k].n1.x, y1 = trg_node[j][k].n1.y, x2 = trg_node[j][k].n2.x, y2 = trg_node[j][k].n2.y;
                            int nstrand = m_msg[x_s].map[y_s].nstrand;
                            if (nstrand == m_msg[x1].map[y1].nstrand ||
                                    m_msg[x_s].map[y_s].nchr != m_msg[x1].map[y1].nchr ||
                                    nstrand * m_msg[x_s].map[y_s].offset < nstrand * m_msg[x2].map[y2].offset ||
                                    nstrand * m_msg[x_e].map[y_e].offset > nstrand * m_msg[x1].map[y1].offset)
                                continue;

                            E_LB(nl, ll) = trg_node[j][k].n1.x;
                            E_RB(nl, ll) = trg_node[j][k].n2.x;
                            L_MF(nl, ll) = L_INTER;
                            if (head == -1) {
                                L_MF(nl, ll) |= L_NMERG;
                                head = l;
                            } else {
                                L_MF(nl, ll) |= L_MERGB;
                                L_MH(nl, ll) = head;
                                _L_MF(line, lsl, head) = L_INTER|L_MERGH;
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
        line_node *mfn = _line_node(line, lsl, m_f[m_ii][0]); int mfl = _line_len(lsl, m_f[m_ii][0]);
        line_node *mfn1 = _line_node(line, lsl, m_f[m_ii+1][0]); int mfl1 = _line_len(lsl, m_f[m_ii+1][0]);
        l = mfn[mfl-1].x - mfn[0].x;
        L = mfn1[mfl1-1].x - mfn1[0].x;
        if (l < 2 && L >= 2) { // DUMP
			for (i = 1; i < m_fn[m_ii]; ++i) {
				_L_MF(line, lsl, m_f[m_ii][i]) = L_DUMP;
                for (_j = ls; _j < ls+len; ++_j) {
                    j = line_rank[_j];
                    line_node *jn = _line_node(line, lsl, j); int jl = _line_len(lsl, j);
                    if (!(L_MF(jn, jl) & L_NMERG) && !(L_MF(jn, jl) & L_MERGH) && !(L_MF(jn, jl) & L_DUMP) && (L_MH(jn, jl) == m_f[m_ii][i]))
                        L_MF(jn, jl) = L_DUMP;
                }
            }
        }
        // right-boundary tri-line
        m_ii = m_i;
        mfn = _line_node(line, lsl, m_f[m_ii][0]); mfl = _line_len(lsl, m_f[m_ii][0]);
        mfn1 = _line_node(line, lsl, m_f[m_ii-1][0]); mfl1 = _line_len(lsl, m_f[m_ii-1][0]);
        l = mfn[mfl-1].x - mfn[0].x;
        L = mfn1[mfl1-1].x - mfn1[0].x;
        if (l < 2 && L >= 2) { // DUMP
            for (i = 1; i < m_fn[m_ii]; ++i) {
                _L_MF(line, lsl, m_f[m_ii][i]) = L_DUMP;
                for (_j = ls; _j < ls+len; ++_j) {
                    j = line_rank[_j];
                    line_node *jn = _line_node(line, lsl, j); int jl= _line_len(lsl, j);
                    if (!(L_MF(jn, jl) & L_NMERG) && !(L_MF(jn, jl) & L_MERGH) && !(L_MF(jn, jl) & L_DUMP) && (L_MH(jn, jl) == m_f[m_ii][i]))
                        L_MF(jn, jl) = L_DUMP;
                }
            }
        }
    }
	for (i = 0; i < len; ++i) { free(m_b[i]); free(m_f[i]); }
    free(m_b); free(m_bn); free(m_f); free(m_fn);
}

void line_filter1(line_node *line, int *lsl, int *line_rank, int ls, int len, int per_max_multi)
{
    int i, _i, j; line_node *nodei; int leni;
	int *m_n = (int*)malloc(len * sizeof(int));
    tri_node **m_b = (tri_node**)malloc(len * sizeof(tri_node*)); // (i, score, NM)
	for (i = 0; i < len; ++i) m_b[i] = (tri_node*)malloc(len * sizeof(tri_node));
    int m_i = -1;
    int b_score, s_score;
    for (_i = ls; _i < ls+len; ++_i) {
        i = line_rank[_i];
        nodei = _line_node(line, lsl, i), leni = _line_len(lsl, i);
		if ((L_MF(nodei, leni) & L_DUMP) || (L_MF(nodei, leni) & L_NMERG)) continue; // dumped and  NOT merged
        else if (L_MF(nodei, leni) & L_MERGH) {    // merged, head
            ++m_i;
            m_b[m_i][0].x = i, m_b[m_i][0].y = L_LS(nodei, leni), m_b[m_i][0].z = L_NM(nodei, leni);
			m_n[m_i] = 1;
        } else { // merged, body
            //if (m_i < 0) { fprintf(stderr, "[line_filter] %s BUG.\n", READ_NAME); exit(1); }
			m_b[m_i][m_n[m_i]].x = i, m_b[m_i][m_n[m_i]].y = L_LS(nodei, leni), m_b[m_i][m_n[m_i]].z = L_NM(nodei, leni);
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
		if (s_score >= b_score/2) f_merge = 1;
		else f_merge = 0;
		// select best & secondary line, set flag
		if (f_merge) {
            // 1. init heap
            ns->node_n = 0;
			for (j = 0; j < m_n[i]; ++j) {
				if (m_b[i][j].y >= b_score/2) {
                    // 2. add node & build heap & update heap
                    int ret = heap_add_node(ns, (line_node){m_b[i][j].x, -1}, m_b[i][j].y, m_b[i][j].z);
                    if (ret == -2) { // dumped
                        _L_MF(line, lsl, m_b[i][j].x) |= L_DUMP;
                    } else if (ret != -1) { // ret is dumped
                        _L_MF(line, lsl, ret) |= L_DUMP;
                    }
                } else _L_MF(line, lsl, m_b[i][j].x) |= L_DUMP;
			} 	
            // 3. ouput_M_nodes, set head/body flags
            build_node_minpos_heap(ns);
            line_node head_n = node_heap_extract_minpos(ns);
            m_head = head_n.x;
            line_node *mhn = _line_node(line, lsl, m_head); int mhl = _line_len(lsl, m_head);
            L_MF(mhn, mhl) = L_MERGH;

            min_l = mhn[0].x;
            max_r = mhn[mhl-1].x;

            int body;
            while ((body = node_heap_extract_minpos(ns).x) != -1) {
                line_node *bdn = _line_node(line, lsl, body); int bdl = _line_len(lsl, body);
                L_MF(bdn, bdl) = L_MERGB;
                L_MH(bdn, bdl) = m_head;
                min_l = MINOFTWO(min_l, bdn[0].x);
                max_r = MAXOFTWO(max_r, bdn[bdl-1].x);
            }
            // set left/right boundary
			L_LB(mhn, mhl) = MINOFTWO(mhn[0].x, min_l);
			L_RB(mhn, mhl) = MAXOFTWO(mhn[mhl-1].x, max_r);
		} else {
			for (j = 0; j < m_n[i]; ++j) {
				if (m_b[i][j].y == b_score)
					_L_MF(line, lsl, m_b[i][j].x) =  L_NMERG;
				else _L_MF(line, lsl, m_b[i][j].x) |= L_DUMP;
			}
		}
	}
    node_free_score(ns);
	for (i = 0; i < len; ++i) free(m_b[i]);
	free(m_b); free(m_n);
}

int line_remove(line_node *line, int *lsl, int *line_rank, int ls, int len)
{
    int l, _l, cur_i=ls;
    // remove line with "DUMP" flag
    for (_l = ls; _l < ls+len; ++_l) {
        l = line_rank[_l];
        line_node *ln = _line_node(line, lsl, l); int ll = _line_len(lsl, l);
        if (!(L_MF(ln, ll) & L_DUMP)) {
            if (cur_i == _l) {
                ++cur_i;
                continue;
            }
            line_rank[cur_i] = l;
            ++cur_i;
        }
    }
    return cur_i-ls;
}

int line_set_bound(line_node *line, int *lsl, int *line_rank, int *line_select_rank,
                   int ls, int *o_len, int left, int right, 
                   trig_node **t_node, int *tri_n, frag_dp_node **f_node, map_msg *m_msg, 
                   int per_max_multi, float ovlp_r) 
{
    if (*o_len <= 0) return 0;
    int i, _i, j, _j, len = *o_len;
    // sort by end-pos
    line_sort_endpos(line, lsl, line_rank, line_select_rank, ls, len);
    // merge lines
    L_MF(line+lsl[line_rank[ls]<<1], lsl[(line_rank[ls]<<1)+1]) = L_NMERG;
    for (i = 1; i < len; ++i)
        line_merge(line_rank[ls+i], line_rank[ls+i-1], line, lsl, ovlp_r);
    // filter best/secondary line, and their tri-nodes, then find inter-line 
    line_filter(line, lsl, line_rank, line_select_rank, ls, len, t_node, tri_n, f_node, m_msg, per_max_multi);
    *o_len = line_remove(line, lsl, line_rank, ls, len);
    len = *o_len;
    // set line-boundary
    int l,r,m,s,e; line_node *nodei; int leni;
    r = right;
    for (_i = ls; _i < ls+len; ++_i) {
        i = line_rank[_i];
        nodei = _line_node(line, lsl, i); 
        leni = _line_len(lsl, i);
        if (!(L_MF(nodei, leni) & L_DUMP)) {
            E_RB(nodei, leni) = r;
            break;
        }
    }
    //if (_i == ls+len) {fprintf(stderr, "[line_set_bound] %s BUG 0.\n", READ_NAME); exit(1);}
    // set right-b
    if (L_MF(nodei, leni) & L_NMERG)
        s = nodei[0].x, e = nodei[leni-1].x;
    else
        s = L_LB(nodei, leni), e = L_RB(nodei, leni);
    m = 0; //Merged head OR NOT-merged head
    for (++_i; _i < ls+len; ++_i) {
        i = line_rank[_i];
        nodei = _line_node(line, lsl, i); leni = _line_len(lsl, i);
        if (L_MF(nodei, leni) & (L_DUMP|L_INTER)) continue;
        if (L_MF(nodei, leni) & L_NMERG) {
            r = s, l = nodei[leni-1].x;
            //set left-b for last-line(s)
            for (_j = _i-1; _j >= m; --_j) {
                j = line_rank[_j];
                line_node *nodej = _line_node(line, lsl, j); int lenj = _line_len(lsl, j);
                if (L_MF(nodej, lenj) & (L_DUMP|L_INTER)) continue;
                E_LB(nodej, lenj) = l;
            }
            s = nodei[0].x, e = nodei[leni-1].x;
            m = _i;
        } else if (L_MF(nodei, leni) & L_MERGH) {
            r = s, l = L_RB(nodei, leni);
            for (_j = _i-1; _j >= m; --_j) {
                j = line_rank[_j];
                if (_L_MF(line, lsl, j) & (L_DUMP|L_INTER)) continue;
                _E_LB(line, lsl, j) = l;
            }
            s = L_LB(nodei, leni), e = L_RB(nodei, leni);
            m = _i;
        } //merged and body
        E_RB(nodei, leni) = r;
    }
    for (_j = _i-1; _j >= m; --_j) {
        j = line_rank[_j];
        if (_L_MF(line, lsl, j) & (L_DUMP|L_INTER)) continue;
        _E_LB(line, lsl, j) = left;
    }
    return 0;
}

int line_set_bound1(line_node *line, int *lsl, int *line_rank, int *line_select_rank,
                   int ls, int *o_len, int left, int right, 
                   int per_max_multi, float ovlp_r) 
{
    if (*o_len <= 0) return 0;
    int i, _i, j, _j, len=*o_len;
    line_node *nodei; int leni;
    // sort by end-pos
    line_sort_endpos(line, lsl, line_rank, line_select_rank, ls, len);
    // merge lines
    _L_MF(line, lsl, line_rank[ls]) = L_NMERG;
    for (i = 1; i < len; ++i)
        line_merge(line_rank[ls+i], line_rank[ls+i-1], line, lsl, ovlp_r);
    // filter best/secondary line, and their tri-nodes, then find inter-line 
    line_filter1(line, lsl, line_rank, ls, len, per_max_multi);
    *o_len = line_remove(line, lsl, line_rank, ls, len);
    len = *o_len;
    // set line-boundary
    int l,r,m,s,e;
    r = right;
    for (_i = ls; _i < ls+len; ++_i) {
        i = line_rank[_i];
        nodei = _line_node(line, lsl, i); leni = _line_len(lsl, i);
        if (!(L_MF(nodei, leni) & L_DUMP)) {
            E_RB(nodei, leni) = r;
            break;
        }
    }
    //if (i == ls+len) {fprintf(stderr, "[line_set_bound] %s BUG 0.\n", READ_NAME); exit(1);}
    // set right-b
    if (L_MF(nodei, leni) & L_NMERG)
        s = nodei[0].x, e = nodei[leni-1].x;
    else
        s = L_LB(nodei, leni), e = L_RB(nodei, leni);
    m = 0; //Merged head OR NOT-merged head
    for (++_i; _i < ls+len; ++_i) {
        i = line_rank[_i];
        nodei = _line_node(line, lsl, i); leni = _line_len(lsl, i);
        if (L_MF(nodei, leni) & (L_DUMP|L_INTER)) continue;
        if (L_MF(nodei, leni) & L_NMERG) {
            r = s, l = nodei[leni-1].x;
            //set left-b for last-line(s)
            for (_j = _i-1; _j >= m; --_j) {
                j = line_rank[_j];
                line_node *nodej = _line_node(line, lsl, j); int lenj = _line_len(lsl, j);
                if (L_MF(nodej, lenj) & (L_DUMP|L_INTER)) continue;
                E_LB(nodej, lenj) = l;
            }
            s = nodei[0].x, e = nodei[leni-1].x;
            m = _i;
        } else if (L_MF(nodei, leni) & L_MERGH) {
            r = s, l = L_RB(nodei, leni);
            for (_j = _i-1; _j >= m; --_j) {
                j = line_rank[_j];
                line_node *nodej = _line_node(line, lsl, j); int lenj = _line_len(lsl, j);
                if (L_MF(nodej, lenj) & (L_DUMP|L_INTER)) continue;
                E_LB(nodej, lenj) = l;
            }
            s = L_LB(nodei, leni), e = L_RB(nodei, leni);
            m = _i;
        } //merged and body
        E_RB(nodei, leni) = r;
    }
    for (_j = _i-1; _j >= m; --_j) {
        j = line_rank[_j];
        line_node *nodej = _line_node(line, lsl, j); int lenj = _line_len(lsl, j);
        if (L_MF(nodej, lenj) & (L_DUMP|L_INTER)) continue;
        E_LB(nodej, lenj) = left;
    }
    return 0;
}

void line_filter_overlap(line_node *line, int *lsl, int *line_rank, int line_n, map_msg *m_msg, frag_dp_node **f_node, int seed_step, int seed_len) 
{
    int i, _i, j, last_i, pre_x, pre_y, cur_x, cur_y;
    line_node *ni; int li;
    for (_i = 0; _i < line_n; ++_i) {
        i = line_rank[_i];
        ni = _line_node(line, lsl, i); li = _line_len(lsl, i);
        last_i = 0;
        for (j = 1; j < li-1; ++j) { // keep the head-tail nodes
            cur_x = ni[j].x, cur_y = ni[j].y;
            pre_x = ni[last_i].x, pre_y = ni[last_i].y;
            //if ((m_msg[cur_x].seed_id - m_msg[pre_x].seed_id) * seed_step < seed_len // read_overlap
            //  || (seed_len + ((m_msg[cur_x].map[cur_y].nstrand==1)?m_msg[pre_x].map[pre_y].len_dif:m_msg[cur_x].map[cur_y].len_dif) > m_msg[cur_x].map[cur_y].nstrand * (m_msg[cur_x].map[cur_y].offset - m_msg[pre_x].map[pre_y].offset) && f_node[cur_x][cur_y].match_flag != F_INSERT)) { // ref_overlap
            if (seed_len + ((m_msg[cur_x].map[cur_y].nstrand==1)?m_msg[pre_x].map[pre_y].len_dif:m_msg[cur_x].map[cur_y].len_dif) > m_msg[cur_x].map[cur_y].nstrand * (m_msg[cur_x].map[cur_y].offset - m_msg[pre_x].map[pre_y].offset) && f_node[cur_x][cur_y].match_flag != F_INSERT) // ref_overlap
                ni[j].x = -1;
            else last_i = j;
        }
        if (li-1 != last_i) {
            cur_x = ni[li-1].x, cur_y = ni[li-1].y;
            pre_x = ni[last_i].x, pre_y = ni[last_i].y;
            //if ((m_msg[cur_x].seed_id - m_msg[pre_x].seed_id) * seed_step < seed_len // read_overlap
            //  || (seed_len + ((m_msg[cur_x].map[cur_y].nstrand==1)?m_msg[pre_x].map[pre_y].len_dif:m_msg[cur_x].map[cur_y].len_dif) > m_msg[cur_x].map[cur_y].nstrand * (m_msg[cur_x].map[cur_y].offset - m_msg[pre_x].map[pre_y].offset) && f_node[cur_x][cur_y].match_flag != F_INSERT)) // ref_overlap
            if (seed_len + ((m_msg[cur_x].map[cur_y].nstrand==1)?m_msg[pre_x].map[pre_y].len_dif:m_msg[cur_x].map[cur_y].len_dif) > m_msg[cur_x].map[cur_y].nstrand * (m_msg[cur_x].map[cur_y].offset - m_msg[pre_x].map[pre_y].offset) && f_node[cur_x][cur_y].match_flag != F_INSERT) // ref_overlap
                ni[last_i].x = -1;
        }
    }
}

void get_fseed_dis(map_msg *m_msg, int pre, int pre_a, int i, int j, int *flag, lamsa_aln_para *AP)
{
    if (pre == -1 || i == -1) { 
        *flag = F_MATCH; 
        return;
    }	//for bound node
    if (pre == i) {
        if (pre_a == j) *flag = F_MATCH; 
        else *flag = F_UNCONNECT;
        return;
    }
    if (m_msg[i].map[j].nchr != m_msg[pre].map[pre_a].nchr || m_msg[i].map[j].nstrand != m_msg[pre].map[pre_a].nstrand)	{ //diff chr or diff strnad
        *flag = F_CHR_DIF; 
        return;
    }

    int seed_len = AP->seed_len; int seed_step = AP->seed_step;
    if (abs(m_msg[pre].seed_id-m_msg[i].seed_id) * seed_step < seed_len) { // for overlapped seeds
        *flag = F_UNCONNECT; return;
    }
    int64_t exp = m_msg[pre].map[pre_a].offset + m_msg[pre].map[pre_a].nstrand * (m_msg[i].seed_id - m_msg[pre].seed_id) * seed_step;	
    int64_t act = m_msg[i].map[j].offset;

    int dis = m_msg[pre].map[pre_a].nstrand * ((m_msg[pre].seed_id < m_msg[i].seed_id)?(act-exp):(exp-act)) - (((m_msg[pre].map[pre_a].nstrand) * (m_msg[pre].seed_id-m_msg[i].seed_id) < 0)?(m_msg[pre].map[pre_a].len_dif):(m_msg[i].map[j].len_dif));

    int mat_dis =  AP->match_dis * ((aln_mode_high_id_err(AP->aln_mode)) ? abs(m_msg[pre].seed_id-m_msg[i].seed_id) : 1); // low-error-rate OR high-error-rate
    if (dis <= mat_dis && dis >= -mat_dis) {
        if (abs(m_msg[pre].seed_id - m_msg[i].seed_id) == 1) *flag = F_MATCH;
        else if (abs(m_msg[pre].seed_id - m_msg[i].seed_id) <= 3*AP->mismatch_thd) *flag = F_MISMATCH; 
        else *flag = F_LONG_MISMATCH;
    } else if (dis > mat_dis && dis < AP->SV_len_thd) *flag = F_DELETE;
    else if ((dis < -mat_dis && dis >= (0-(abs(m_msg[i].seed_id-m_msg[pre].seed_id)*seed_step-seed_len))) // nonoverlaped ins
          || (dis < -(AP->split_len/2) && dis >= -AP->SV_len_thd)) { // overlaped ins/dup
        *flag = F_INSERT; 
    } else { 
        *flag = F_UNCONNECT; 
        return;
    }
}

void fnode_set(frag_dp_node *f_node, line_node from, 
               int score, int NM,
               uint8_t match_flag, int dp_flag)
{
    (*f_node).son_flag = F_INIT;
    (*f_node).from = from;
    (*f_node).score = score;
    (*f_node).tol_NM = NM;

    (*f_node).match_flag = match_flag;
    (*f_node).dp_flag = dp_flag;

    // general init
    (*f_node).node_n = 1;
    (*f_node).in_de = 0;
    (*f_node).son_n = 0;
    (*f_node).max_score = score;
    (*f_node).max_NM = NM;
    (*f_node).max_node.x = (*f_node).seed_i, (*f_node).max_node.y = (*f_node).aln_i;
}

int frag_dp_init(frag_dp_node **f_node, 
                 map_msg *m_msg, int seed_i, 
                 line_node from, 
                 lamsa_aln_para *AP,
                 int dp_flag)
{
    int i;
    if (from.x == START_NODE.x) { //UNLIMITED
        for (i = 0; i < m_msg[seed_i].map_n; ++i)
            //init-score: 0 == 1 + 1 - F_SV_PEN
            fnode_set(&(f_node[seed_i][i]), from, 1, m_msg[seed_i].map[i].NM, F_MATCH, dp_flag);
    } else {
        int con_flag, con_score;
        for (i = 0; i < m_msg[seed_i].map_n; ++i) {
            
            get_fseed_dis(m_msg, from.x, from.y, seed_i, i, &con_flag, AP);
            con_score = AP->frag_score_table[con_flag];

            if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
                fnode_set(&(f_node[seed_i][i]), from, 2+con_score, m_msg[from.x].map[from.y].NM+m_msg[seed_i].map[i].NM, con_flag, dp_flag);
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
int frag_dp_update(frag_dp_node **f_node, map_msg *m_msg, 
                   int seed_i, int aln_i, int start, 
                   lamsa_aln_para *AP, int dp_flag)
{
    int i, j, con_flag, con_score;
    line_node max_from;
    int max_score, max_flag, max_NM;

    max_from = f_node[seed_i][aln_i].from;
    max_score = f_node[seed_i][aln_i].score;
    max_NM = f_node[seed_i][aln_i].tol_NM; max_flag = f_node[seed_i][aln_i].dp_flag;

    for (i = seed_i - 1; i >= start; --i) {
        for (j = 0; j < m_msg[i].map_n; ++j) {
            if (f_node[i][j].dp_flag == dp_flag)
            {
                // '+' strand: successor-Match
                if (m_msg[i].map[j].nstrand == 1) { 
                    if (f_node[i][j].son_flag <= F_MATCH_THD) continue;
                }

                get_fseed_dis(m_msg, i, j, seed_i, aln_i, &con_flag, AP);
                if (con_flag == F_UNCONNECT || con_flag == F_CHR_DIF) continue;
                
                // '-' strand: precursor-Match 
                if (m_msg[i].map[j].nstrand == -1 && con_flag <= F_MATCH_THD) { 
                    max_from = (line_node){i,j};
                    con_score = AP->frag_score_table[con_flag];
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_NM = f_node[i][j].tol_NM + f_node[seed_i][aln_i].tol_NM;
                    goto UPDATE;
                }

                con_score = AP->frag_score_table[con_flag];
                if (f_node[i][j].score + 1 + con_score > max_score)	{
                    max_from = (line_node){i,j};
                    max_score = f_node[i][j].score + 1 + con_score;
                    max_flag = con_flag;
                    max_NM = f_node[i][j].tol_NM + f_node[seed_i][aln_i].tol_NM;
                } else if (f_node[i][j].score + 1 + con_score == max_score) {
                    if (f_node[seed_i][aln_i].tol_NM + f_node[i][j].tol_NM < max_NM) {
                        max_from = (line_node){i,j};
                        max_score = f_node[i][j].score + 1 + con_score;
                        max_flag = con_flag;
                        max_NM = f_node[i][j].tol_NM + f_node[seed_i][aln_i].tol_NM;
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
        f_node[seed_i][aln_i].tol_NM = max_NM;
        f_node[seed_i][aln_i].match_flag = max_flag;
        f_node[seed_i][aln_i].node_n = f_node[max_from.x][max_from.y].node_n+1;

        fnode_add_son(f_node, max_from, (line_node){seed_i, aln_i}); //generate tree
    }
    return 0;
}

int frag_dp_per_init(frag_dp_node **f_node, map_msg *m_msg, 
                     int seed_i, int aln_i, line_node from, 
                     lamsa_aln_para *AP,
                     int dp_flag)
{
    if (from.x == START_NODE.x) //UNLIMITED
        //init-score: 0 == 1 + 1 - F_SV_PEN
        fnode_set(&(f_node[seed_i][aln_i]), from, 1, m_msg[seed_i].map[aln_i].NM, F_MATCH, dp_flag);
    else {
        int con_flag, con_score;
        get_fseed_dis(m_msg, from.x, from.y, seed_i, aln_i, &con_flag, AP);
        con_score = AP->frag_score_table[con_flag];

        if (con_flag != F_UNCONNECT && con_flag != F_CHR_DIF)
            fnode_set(&(f_node[seed_i][aln_i]), from, 2+con_score, m_msg[seed_i].map[aln_i].NM+m_msg[from.x].map[from.y].NM, con_flag, dp_flag);
        else f_node[seed_i][aln_i].dp_flag = 0-dp_flag;
    }
    return 0;
}

void node_add_score(int score, int NM, line_node node, node_score *ns, frag_dp_node **f_node)
{
    if (score < ns->min_score_thd) return;
    if (ns->node_n <= ns->max_n-1) {
        ns->score[ns->node_n] = score;
        ns->NM[ns->node_n] = NM;
        ns->node[(ns->node_n)++] = node;
        f_node[node.x][node.y].dp_flag = TRACKED_FLAG;
        line_node tmp = f_node[node.x][node.y].from;
        while (tmp.x != -1) {
            f_node[tmp.x][tmp.y].dp_flag = TRACKED_FLAG;
            tmp = f_node[tmp.x][tmp.y].from;
        }
    } else {
        fprintf(stderr, "[lamsa_aln] node_add_score ERROR. (%d %d)\n", ns->node_n, ns->max_n); exit(1);
    }
}

//1.F_MATCH/F_MISMATCH/F_LONG_MISMATCH
//2.F_INS/F_DEL:
//  a.Better-Score
//  b.Short-Distance
line_node get_max_son(frag_dp_node **f_node, int x, int y)
{
    int i, max_score, max_NM, max_flag, max_dis, flag_thd=F_INIT;
    line_node max, son;
    max_score = 0, max_NM = 0, max_dis = 0;

    for (i = 0; i < f_node[x][y].son_n; ++i) {
        son = f_node[x][y].son[i];
        max_flag = f_node[son.x][son.y].match_flag;
        // for high-err-overlap seed
        //if (max_flag == F_MATCH || max_flag == F_MISMATCH)// || max_flag == F_LONG_MISMATCH)
        //    return son;
        if (f_node[son.x][son.y].match_flag <= flag_thd && (f_node[son.x][son.y].max_score > max_score || (f_node[son.x][son.y].max_score == max_score && (son.x - x < max_dis || f_node[son.x][son.y].max_NM < max_NM)))) {
            max = son;
            max_score = f_node[son.x][son.y].max_score;
            max_NM = f_node[son.x][son.y].max_NM;
            max_dis = son.x - x;
            if (f_node[son.x][son.y].match_flag <= F_MATCH_THD) flag_thd = F_MATCH_THD;
        } 
    }
    return max;
}

void cut_branch(frag_dp_node **f_node, map_msg *m_msg, int x, int y, node_score *ns) {
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
        f_node[son.x][son.y].max_NM -= (f_node[son.x][son.y].tol_NM - m_msg[son.x].map[son.y].NM);
        max_node = f_node[son.x][son.y].max_node;
        f_node[max_node.x][max_node.y].node_n -= (f_node[son.x][son.y].node_n-1);
        node_add_score(f_node[son.x][son.y].max_score, f_node[son.x][son.y].max_NM, max_node , ns, f_node);
    }
    // refresh max_socre/node of (x,y), based on max_son, consider negative edge
    if (f_node[x][y].score > f_node[max_son.x][max_son.y].max_score) { // negative edge
        f_node[max_son.x][max_son.y].in_de = -1;
        f_node[max_son.x][max_son.y].from = START_NODE;
        f_node[max_son.x][max_son.y].max_score -= (f_node[max_son.x][max_son.y].score - 1);
        f_node[max_son.x][max_son.y].max_NM -= (f_node[max_son.x][max_son.y].tol_NM - m_msg[max_son.x].map[max_son.y].NM);
        max_node = f_node[max_son.x][max_son.y].max_node;
        f_node[max_node.x][max_node.y].node_n -= (f_node[max_son.x][max_son.y].node_n-1);
        node_add_score(f_node[max_son.x][max_son.y].max_score, f_node[max_son.x][max_son.y].max_NM, max_node, ns, f_node);
        f_node[x][y].son_n = 0;
        f_node[x][y].max_node = (line_node){x,y};
        f_node[x][y].max_score = f_node[x][y].score;
        f_node[x][y].max_NM = f_node[x][y].tol_NM;
    } else {
        f_node[x][y].son_n = 1;
        f_node[x][y].son[0] = max_son;
        f_node[x][y].max_node = f_node[max_son.x][max_son.y].max_node;
        f_node[x][y].max_score = f_node[max_son.x][max_son.y].max_score;
        f_node[x][y].max_NM = f_node[max_son.x][max_son.y].max_NM;
    }
    f_node[x][y].in_de = 0; // leaf-node, avaible for branch-track
}

//backtrack from (x,y)
void branch_track_new(frag_dp_node **f_node, map_msg *m_msg, int x, int y, node_score *ns)
{
    int max_score, max_NM; line_node max_node; 
    f_node[x][y].in_de = -1;
    if (f_node[x][y].son_n == 0) {
        max_node = f_node[x][y].max_node = (line_node){x,y};
        max_score = f_node[x][y].max_score = f_node[x][y].score;
        max_NM = f_node[x][y].max_NM = f_node[x][y].tol_NM;
    } else {
        max_node = f_node[x][y].max_node;
        max_score = f_node[x][y].max_score;
        max_NM = f_node[x][y].max_NM;
    }

    line_node fa = f_node[x][y].from, son; int fa_son_n;

    while (fa.x != START_NODE.x) {
        fa_son_n = f_node[fa.x][fa.y].son_n;
        if (fa_son_n == 1) { // backtrack
            if (f_node[fa.x][fa.y].score > max_score) { // negative edge
                son = f_node[fa.x][fa.y].son[0];
                f_node[son.x][son.y].in_de = -1;
                f_node[son.x][son.y].from = START_NODE;
                f_node[son.x][son.y].max_score -= (f_node[son.x][son.y].score - 1);
                f_node[son.x][son.y].max_NM -= (f_node[son.x][son.y].tol_NM - m_msg[son.x].map[son.y].NM);
                f_node[max_node.x][max_node.y].node_n -= (f_node[son.x][son.y].node_n-1);
                node_add_score(f_node[son.x][son.y].max_score, f_node[son.x][son.y].max_NM, max_node, ns, f_node);
                f_node[fa.x][fa.y].son_n = 0;
                max_score = f_node[fa.x][fa.y].score;
                max_NM = f_node[fa.x][fa.y].tol_NM;
                max_node = fa;
            }
            f_node[fa.x][fa.y].max_score = max_score;
            f_node[fa.x][fa.y].max_NM = max_NM;
            f_node[fa.x][fa.y].max_node = max_node;
            f_node[fa.x][fa.y].in_de = -1;
            fa = f_node[fa.x][fa.y].from;
        } else { // fa_son_n > 1, TREE_BRANCH
            --f_node[fa.x][fa.y].in_de;
            if (f_node[fa.x][fa.y].in_de == 0)
                cut_branch(f_node, m_msg, fa.x, fa.y, ns);
            return;
        }
    }
    // TREE_START; 
    node_add_score(max_score, max_NM, max_node, ns, f_node);// add (max_node, max_score);
    return;
}

//for multi-dp-lines
int frag_mini_dp_multi_line(frag_dp_node **f_node, map_msg *m_msg, 
                            lamsa_aln_para *AP, kseq_t *seqs,
                            int left_b, int right_b, 
                            reg_t trg_reg,
                            line_node *line, int *lsl,
							int line_n_max)
{
    if (left_b+1 >= right_b) return 0;

    line_node head = START_NODE; // unlimited
    int start, end;
    int i, j, dp_flag = WHOLE_FLAG;// (whole DP)
    int l_i, node_i;
    line_node left, right;
    line_node _right;
    line_node *node; int len;
    //int candi_n; line_node candi[max_multi];

    left.x = left_b, right.x = right_b;
    start = left.x+1;
    end = right.x-1;

    //first dp int
    for (i = start; i <= end; ++i) {
        for (j = 0; j < m_msg[i].map_n; ++j) {
            if (f_node[i][j].dp_flag != TRACKED_FLAG)
                frag_dp_per_init(f_node, m_msg, i, j, head, AP, dp_flag);
        }
    }
    //dp update
    for (i = start+1; i <= end; ++i) {
        for (j = 0; j < m_msg[i].map_n; ++j) {
            if (f_node[i][j].dp_flag == dp_flag)
                frag_dp_update(f_node, m_msg, i, j, start, AP, dp_flag);
        }
    }
    // tree-pruning
    node_score *ns = node_init_score(line_n_max); //size of node_score
    ns->min_score_thd = 0;
    //node_score *ns = node_init_score(80); //size of node_score
    for (i = end; i >=start; --i) {
        for (j = 0; j < m_msg[i].map_n; ++j) {
            if (f_node[i][j].dp_flag == dp_flag && f_node[i][j].in_de == 0) //leaf node
                branch_track_new(f_node, m_msg, i, j, ns); 
        }
    }
    l_i = 0;
    int score, NM;
    int next_start = 0;
    while (1) {
        _right = node_pop(ns, &score, &NM);
        if (_right.x == START_NODE.x) break;
        node_i = f_node[_right.x][_right.y].node_n-1;
        lsl[l_i<<1] = next_start; lsl[(l_i<<1)+1] = node_i+1;
        next_start += (node_i+1+L_EXTRA);
        int s_i = f_node[_right.x][_right.y].seed_i, a_i = f_node[_right.x][_right.y].aln_i;
        for (i = 0; i < trg_reg.beg_n; i++) {
            if (m_msg[s_i].map[a_i].nchr == trg_reg.ref_beg[i].chr
                && labs((m_msg[s_i].map[a_i].offset-trg_reg.ref_beg[i].ref_pos) - (s_i-left_b)*AP->seed_step) < AP->SV_len_thd) {
                // connectable with main-line
                if (score > 1) score += score/2;
                else score++;
                goto L_START;
            }
        }
        for (i = 0; i < trg_reg.end_n; i++) {
            if (m_msg[s_i].map[a_i].nchr == trg_reg.ref_end[i].chr
                && labs((m_msg[s_i].map[a_i].offset-trg_reg.ref_end[i].ref_pos) - (s_i-left_b)*AP->seed_step) < AP->SV_len_thd) {
                // connectable with main-line
                if (score > 1) score += score/2;
                else score++;
                goto L_START;
            }
        }
L_START:
        node = _line_node(line, lsl, l_i); 
        len = _line_len(lsl, l_i);
        L_LS(node, len) = score;
        L_BS(node, len) = score;
        L_NM(node, len) = NM;
        while (_right.x != head.x) {
            if (node_i < 0) { 
                fprintf(stderr, "\n[frag mini dp multi] %s node_i BUG 1.\n", seqs->name.s); exit(1); 
			}
            node[node_i--] = _right;
            _right = f_node[_right.x][_right.y].from;
        }
        if (node_i >= 0) { 
            fprintf(stderr, "\n[frag mini dp multi] %s node_i BUG 2.\n", seqs->name.s); exit(1); 
		}
        ++l_i;
    }
    node_free_score(ns);
    return l_i;
}

int trg_dp_line(frag_dp_node **f_node, map_msg *m_msg,
                lamsa_aln_para *AP, kseq_t *seqs,
                int left, int right, 
                reg_t trg_reg,
                line_node *line, int *lsl, int *line_rank, int *line_select_rank,
                int line_n_max) 
{
    int l = frag_mini_dp_multi_line(f_node, m_msg, AP, seqs, left, right, trg_reg, line, lsl, line_n_max);//, 0, 0);
    line_set_bound1(line, lsl, line_rank, line_select_rank, 0,  &l, left, right, AP->ske_max, AP->ovlp_rat);
    return l;
}

void frag_min_extend(frag_dp_node **f_node, map_msg *m_msg,
                    int node_i, int aln_i,
                    int aln_min, int dp_flag,
                    lamsa_aln_per_para *APP, lamsa_aln_para *AP)
{
    int i, j, con_flag;
    int last_x = node_i, last_y = aln_i;
    // 0 -> node_i-1, node_i+1 -> node_n-1
    //from right to left 
    i = last_x - 1;
    while (i >= 0) {
        if (m_msg[i].map_n > aln_min) {
            for (j = 0; j < m_msg[i].map_n; ++j) {
                get_fseed_dis(m_msg, i, j, last_x, last_y, &con_flag, AP);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH || con_flag == F_LONG_MISMATCH) {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        --i;
    }
    i = last_x + 1;
    while (i < APP->seed_out) {
        if (m_msg[i].map_n > aln_min) {
            for (j = 0; j < m_msg[i].map_n; ++j) {
                get_fseed_dis(m_msg, last_x, last_y, i, j, &con_flag, AP);
                if (con_flag == F_MATCH || con_flag == F_MISMATCH || con_flag == F_LONG_MISMATCH) {
                    f_node[i][j].dp_flag = dp_flag;
                    break;
                }
            }
        }
        ++i;
    }
}

int frag_mini_dp_line(frag_dp_node **f_node, map_msg *m_msg, 
                      lamsa_aln_para *AP, kseq_t *seqs,
                      line_node left, line_node right, 
                      line_node *line, int *de_score, int *de_NM,
                      int _head, int _tail)
{
    extern int f_BCC_score_table[10];
    line_node head = ((_head)?left:START_NODE);
    int old_score, old_NM;
    if (_tail==0) {
        old_score = 1;
        old_NM = (left.x==START_NODE.x)?0:m_msg[left.x].map[left.y].NM;
    } else {
        old_score = 2+f_BCC_score_table[f_node[right.x][right.y].match_flag];
        old_NM = ((left.x==START_NODE.x)?0:m_msg[left.x].map[left.y].NM) + m_msg[right.x].map[right.y].NM;
    }
    int i, j, dp_flag = MULTI_FLAG;
    //dp init
    for (i = left.x + 1; i < right.x; ++i) {
        for (j = 0; j < m_msg[i].map_n; ++j) {
            if (f_node[i][j].dp_flag == dp_flag || f_node[i][j].dp_flag == 0-dp_flag)
                frag_dp_per_init(f_node, m_msg, i, j, head, AP, dp_flag);
        }
    }
    //dp update
    for (i = left.x + 2; i < right.x; ++i) {
        for (j = 0; j < m_msg[i].map_n; ++j) {
            if (f_node[i][j].dp_flag == dp_flag)
                frag_dp_update(f_node, m_msg, i, j, left.x+1, AP, dp_flag);
        }
    }
    //find backtrack start node
    int max_score, max_NM = 0, max_n = 0;
    line_node max_node = head;//left;
    //UNLIMITED 
    if (_tail == 0) {
        max_score = old_score;
        for (i = right.x - 1; i > left.x; --i) {
            for (j = 0; j < m_msg[i].map_n; ++j) {
                if (f_node[i][j].dp_flag == dp_flag) {
                    if (f_node[i][j].score > max_score) {
                        max_score = f_node[i][j].score;
                        max_NM = f_node[i][j].tol_NM;
                        max_node = (line_node){i, j};
                        max_n = f_node[i][j].node_n;
                    } else if (f_node[i][j].score == max_score) {
                        if (f_node[i][j].tol_NM < max_NM) {
                            max_score = f_node[i][j].score;
                            max_NM = f_node[i][j].tol_NM;
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
        f_node[right.x][right.y].tol_NM = old_NM;
        f_node[right.x][right.y].node_n = 1;
        frag_dp_update(f_node, m_msg, right.x, right.y, left.x+1, AP, dp_flag);
        max_score = f_node[right.x][right.y].score;
        max_NM = f_node[right.x][right.y].tol_NM;
        max_node = f_node[right.x][right.y].from;
        max_n = f_node[right.x][right.y].node_n-1; // right is NOT inculding 
    }
    //backtrack
    line_node _right = max_node;
    int node_i = max_n-1;
    //while (_right.x != left.x)
    while (_right.x != head.x) {
        if (node_i < 0) { 
            fprintf(stderr, "\n[frag mini dp] %s node_i BUG 1.\n", seqs->name.s); exit(1); 
        }
        line[node_i--] = _right;
        _right = f_node[_right.x][_right.y].from;
    }
    if (node_i >= 0) { fprintf(stderr, "\n[frag mini dp] %s node_i BUG 2.\n", seqs->name.s); exit(1); }
    *de_score += (max_score-old_score);
    *de_NM += (max_NM-old_NM);
    return max_n;
}

int frag_dp_path(map_msg *m_msg, frag_msg **f_msg,
                 lamsa_aln_para *AP, lamsa_aln_per_para *APP, 
                 int *l_n, int *line_m,
                 line_node *line, int *lsl, int *line_rank,
                 frag_dp_node ***f_node)
{
    if (*l_n == 0) return 0;
    int line_n = *l_n, i, l, _l;

#ifdef __DEBUG__
    int _i;
    for (_i = 0; _i < line_n; ++_i) {
        i = line_rank[_i];
        line_node *ni = _line_node(line, lsl, i);
        fprintf(stderr, "%d(%d,%d score: %d):\t", i+1, _L_MF(line, lsl, i), _L_MH(line, lsl, i), _L_LS(line, lsl, i));
        for (l = 0; l < _line_len(lsl, i); ++l)
            if (ni[l].x >=0 ) fprintf(stderr, "(%d, %d)\t", m_msg[ni[l].x].seed_id-1, ni[l].y);
        fprintf(stderr, "\n");
    }
#endif

    int frag_num;
    int cur_x , cur_y , pre_x, pre_y;

    int left_bound, right_bound;

    if (line_n > (*line_m)) {
        (*f_msg) = (frag_msg*)realloc(*f_msg, line_n * sizeof(frag_msg));
        for (i = (*line_m); i < line_n; ++i) {
            frag_init_msg((*f_msg)+i, (*f_msg)->frag_max); 
            frag_copy_msg(*f_msg, (*f_msg)+i);
        }
        if ((*f_msg) == NULL) { fprintf(stderr, "\n[frag_dp_path] Not enough memory.(line_m: %d)\n", line_n); exit(1); }
        (*line_m)= line_n;
    }

    if (aln_mode_overlap_seed(AP->aln_mode)) line_filter_overlap(line, lsl, line_rank, line_n, m_msg, *f_node, AP->seed_step, AP->seed_len);

    for (_l = 0; _l < line_n; ++_l) {
        l = line_rank[_l];
        line_node *ln = _line_node(line, lsl, l); int ll = _line_len(lsl, l);
        frag_num=0;
        pre_x = ln[ll-1].x; pre_y = ln[ll-1].y;
        right_bound = APP->seed_all+1;
        frag_set_msg(m_msg, pre_x, pre_y, FRAG_END, (*f_msg)+_l, frag_num);
        ((*f_msg)+_l)->frag_right_bound = right_bound;
        for (i = ll-1; i > 0; --i) {
            cur_x = pre_x; cur_y = pre_y;
            if (ln[i-1].x < 0) continue;
            pre_x = ln[i-1].x; pre_y = ln[i-1].y;
            //INS
            if ((*f_node)[cur_x][cur_y].match_flag == F_INSERT)
            {
                frag_set_msg(m_msg, cur_x, cur_y, FRAG_START, (*f_msg)+_l, frag_num);
                ++frag_num;
                frag_set_msg(m_msg, pre_x, pre_y, FRAG_END, (*f_msg)+_l, frag_num);
            }
            //DEL
            else if ((*f_node)[cur_x][cur_y].match_flag == F_DELETE) {
                frag_set_msg(m_msg, cur_x, cur_y, FRAG_START, (*f_msg)+_l, frag_num);
                ++frag_num;
                frag_set_msg(m_msg, pre_x, pre_y, FRAG_END, (*f_msg)+_l, frag_num);
            }
            //MIS	
            else if ((*f_node)[cur_x][cur_y].match_flag == F_MISMATCH || (*f_node)[cur_x][cur_y].match_flag == F_LONG_MISMATCH) {
                //take mis-match as NEW-FLAG case
                //frag_set_msg(m_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+cur_line, frag_num, seed_len);
                frag_set_msg(m_msg, cur_x, cur_y, FRAG_START, (*f_msg)+_l, frag_num);
                ++frag_num;
                frag_set_msg(m_msg, pre_x, pre_y, FRAG_END, (*f_msg)+_l, frag_num);
                //find new flag in mis-match seeds
                //dp-line for INV/TRS
            } else if ((*f_node)[cur_x][cur_y].match_flag == F_MATCH)	//: MATCH
            {
                frag_set_msg(m_msg, pre_x, pre_y, FRAG_SEED, (*f_msg)+_l, frag_num);
            } else { fprintf(stderr, "\n[frag dp path] Error: Unknown flag, \"%d\"\n", (*f_node)[cur_x][cur_y].match_flag); exit(1);}
        }
        //last start
        cur_x = ln[0].x; cur_y = ln[0].y;
        frag_set_msg(m_msg, cur_x, cur_y, FRAG_START, (*f_msg)+_l, frag_num);
        left_bound = 0;

        //MIS-MATCH
        ((*f_msg)+_l)->frag_left_bound = left_bound;
        ((*f_msg)+_l)->line_score = L_LS(ln, ll);
    }
    *l_n = _l;

#ifdef __DEBUG__
    fprintf(stderr, "After Filtering:\n");
    for (_i = 0; _i < line_n; ++_i) {
        i = line_rank[_i];
        line_node *ni = _line_node(line, lsl, i);
        fprintf(stderr, "%d(%d,%d score: %d):\t", i+1, _L_MF(line, lsl, i), _L_MH(line, lsl, i), _L_LS(line, lsl, i));
        for (l = 0; l < _line_len(lsl, i); ++l) {
            if (ni[l].x >=0 ) fprintf(stderr, "(%d, %d)\t", m_msg[ni[l].x].seed_id-1, ni[l].y);
        }
        fprintf(stderr, "\n");
    }
#endif

    return 1;
}

int frag_line_remain(aln_reg *a_reg, map_msg *m_msg, frag_msg **f_msg,
                     lamsa_aln_per_para *APP, lamsa_aln_para *AP, kseq_t *seqs,
                     line_node *line, int *lsl, int *line_rank, int *line_select_rank,
                     int *line_m, frag_dp_node ***f_node, 
                     line_node *_line, int *_lsl, int *_line_rank,
                     int line_n_max)
{
    int l_n = 0;
    aln_reg *re_reg = aln_init_reg(seqs->seq.l);
    if (get_remain_reg(a_reg, re_reg, AP, AP->seed_len, seqs->seq.l) == 0) goto End;
    int i, j, _j, k;
    int left_id, right_id, left, right, l;
    int next_start = 0;
    for (i = 0; i < re_reg->reg_n; ++i) {
        // get region boundary (left, right)
        left_id = (re_reg->reg[i].beg + AP->seed_inv - 1) / AP->seed_step+1; 
        right_id = (re_reg->reg[i].end - 1) / AP->seed_step+1;
        if (right_id > APP->seed_all) right_id -= 1;
        left = right = -2;
        for (j = 0; j < APP->seed_out; ++j) {
            if (m_msg[j].seed_id >= left_id) {
                left = j-1;
                break;
            }
        }
        if (left == -2) continue;
        for (j = APP->seed_out-1; j >= 0; --j) {
            if (m_msg[j].seed_id <= right_id) {
                right = j+1;
                break;
            }
        }
        if (right == -2) continue;
        // store line-info in _line
        l = trg_dp_line(*f_node, m_msg, AP, seqs, left, right, re_reg->reg[i], _line, _lsl, _line_rank, line_select_rank, line_n_max);
        // copy _line, _lsl, _line_rank
        for (_j = 0; _j < l; ++_j) {
            j = _line_rank[_j];
            lsl[(l_n+_j)<<1] = next_start; lsl[((l_n+_j)<<1)+1] = _lsl[(j<<1)+1];
            for (k = 0; k < _lsl[(j<<1)+1]+L_EXTRA; ++k)
                line[next_start+k] = _line[_lsl[(j<<1)]+k];
            next_start += (_lsl[(j<<1)+1] + L_EXTRA);
            line_rank[l_n+_j] = l_n+_j;
        }
        l_n += l;
    }
End:
    aln_free_reg(re_reg);
    frag_dp_path(m_msg, f_msg, AP, APP, &l_n, line_m, line, lsl, line_rank, f_node);
    return l_n;
}

//new tree-generating and pruning
int frag_line_BCC(map_msg *m_msg, frag_msg **f_msg,
        lamsa_aln_per_para *APP, lamsa_aln_para *AP, kseq_t *seqs,
        line_node *line, int *lsl, int *line_rank, int *line_select_rank, int *line_m,
        frag_dp_node ***f_node,
        line_node *_line, int line_n_max)
{
    int i, j, k, line_n; line_node *ln;
    int min_n = AP->first_loci_thd, min_exist=0, min_num = 0;

    { // dp init
        for (i = 0; i < APP->seed_out; ++i)
        {
            if (m_msg[i].map_n <= min_n)
            {
                frag_dp_init(*f_node, m_msg, i, START_NODE, AP, MIN_FLAG);
                min_exist = 1;
                ++min_num;
            } else frag_dp_init(*f_node, m_msg, i, START_NODE, AP, MULTI_FLAG);
        }
        if (!min_exist || min_num * 3 < APP->seed_out)
        {
            for (i = 0; i < APP->seed_out; ++i) {
                for (j = 0; j < m_msg[i].map_n; ++j) (*f_node)[i][j].dp_flag = MIN_FLAG;
            }
            min_n = AP->per_aln_m;
            min_exist = 1;
        }
    }
    { // dp update
        //min extend, when min_n == PER_ALN_N: no need to extend
        if (min_n != AP->per_aln_m)
        {
            for (i = 0; i < APP->seed_out; ++i) {
                if (m_msg[i].map_n <= min_n) {
                    for (j = 0; j < m_msg[i].map_n; ++j)
                        frag_min_extend(*f_node, m_msg, i, j, min_n, MIN_FLAG, APP, AP);
                }
            }
        }
        //min update
        for (i = 1; i < APP->seed_out; ++i) {
            for (j = 0; j < m_msg[i].map_n; ++j) {
                if ((*f_node)[i][j].dp_flag == MIN_FLAG)
                    frag_dp_update(*f_node, m_msg, i, j, 0/*update start pos*/, AP, MIN_FLAG);
            }
        }
    }
    { // backtrack
        //prune tree to find multi-nodes/roads for backtrack
        node_score *ns = node_init_score(line_n_max); //size of node_score
        ns->min_score_thd = 2;
        for (i = APP->seed_out-1; i >=0; --i) {
            for (j = 0; j < m_msg[i].map_n; ++j) {
                if ((*f_node)[i][j].dp_flag == MIN_FLAG && (*f_node)[i][j].in_de == 0) //leaf node
                    branch_track_new(*f_node, m_msg, i, j, ns); 
            }
        }
        int line_score, line_NM; line_node max_node;
        int l_i = 0, min_l=ns->node_n, o_l = min_l;
        int node_i, mini_len, next_start=0;
        trig_node **inter_tri = (trig_node**)malloc(o_l* sizeof(trig_node*)); int *inter_n=(int*)malloc(o_l*sizeof(int));
        for (i = 0; i < o_l;++i)
            inter_tri[i] = (trig_node*)malloc(APP->seed_out * sizeof(trig_node));
        //multi-backtrack
        //extend for min-line AND set inter-trigger
        while (1) {
            //pop max score in the heap
            max_node = node_pop(ns, &line_score, &line_NM);    
            if (max_node.x == -1) break;
            //backtrack for max-node
            line_node last_n, right, left;
            //fill-up for right end
            node_i=0;
            inter_n[l_i]=0, lsl[l_i<<1] = next_start, ln = line+next_start;
            if (max_node.x < APP->seed_out-1) {
                mini_len = frag_mini_dp_line(*f_node, m_msg, AP, seqs, max_node, (line_node){APP->seed_out, 0}, _line, &line_score, &line_NM, 1, 0);
                for (k = mini_len-1; k >= 0; --k) {
                    ln[node_i++] = _line[k];
                    (*f_node)[_line[k].x][_line[k].y].dp_flag = TRACKED_FLAG; 
                }
                ln[node_i] = max_node;
                //set inter-trigger
                last_n = ln[0];
                for (k = mini_len-1; k >= 0; --k) {
                    if (last_n.x - ln[node_i-k].x > 2) {
                        // pre-trigger: (line[l_i][node_i-k].x , right.x)
                        inter_tri[l_i][inter_n[l_i]++] = (trig_node){(line_node){ln[node_i-k].x, ln[node_i-k].y}, (line_node){last_n.x, last_n.y}, 2};
                    }
                    last_n = ln[node_i-k];
                }
            }
            right = max_node;
            // fill-up for every gap
            while (right.x != START_NODE.x) {
                ln[node_i++] = right;
                left = (*f_node)[right.x][right.y].from;
                if (left.x < right.x - 1 )
                {
                    //mini-dp with anchors
                    mini_len = frag_mini_dp_line(*f_node, m_msg, AP, seqs, left, right, _line, &line_score, &line_NM, 1, 1);
                    for (k = mini_len-1; k >= 0; --k) {
                        ln[node_i++] = _line[k];
                        (*f_node)[_line[k].x][_line[k].y].dp_flag = TRACKED_FLAG; 
                    }
                    ln[node_i] = left;	//twice write, for the last multi-dp-line
                    //set inter-trigger
                    last_n = right;
                    for (k = mini_len; k >= 0; --k) {
                        if (last_n.x - ln[node_i-k].x > 2) {
                            // pre-trigger: ((line+line_start[l_i][node_i-k].x , right.x)
                            if (ln[node_i-k].x == START_NODE.x) continue;
                            inter_tri[l_i][inter_n[l_i]++] = (trig_node){(line_node){ln[node_i-k].x, ln[node_i-k].y},(line_node){last_n.x, last_n.y}, 2};
                        }
                        last_n = ln[node_i-k];
                    }
                }
                right = left;
            }
            //invert line
            line_node tmp;
            for (k = 0; k < node_i/2; ++k) {
                tmp =ln[k]; ln[k] = ln[node_i-k-1]; ln[node_i-k-1] = tmp;
            }
            lsl[(l_i<<1)+1] = node_i;

            L_LS(ln, node_i) = line_score; L_BS(ln, node_i) = line_score; L_NM(ln, node_i) = line_NM;
            l_i++; next_start += (node_i+L_EXTRA);
        }
        // min-lines have been extended, and sorted by end-pos
        // set boundary and filter best/secondary
        line_set_bound(line, lsl, line_rank, line_select_rank, 0, &min_l, -1, APP->seed_out, inter_tri, inter_n, *f_node, m_msg, AP->ske_max, AP->ovlp_rat);
        for (i = 0; i < o_l; ++i) free(inter_tri[i]);
        free(inter_tri); free(inter_n);
        node_free_score(ns);

        line_n = min_l;

        frag_dp_path(m_msg, f_msg, AP, APP, &line_n, line_m, line, lsl, line_rank, f_node);
        return line_n;
    }
}
