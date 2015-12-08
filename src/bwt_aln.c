#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "frag_check.h"
#include "lamsa_aln.h"
#include "split_mapping.h"
#include "bwt_aln.h"
#include "lamsa_heap.h"
#include "bntseq.h"
#include "bwt.h"
#include "ksw.h"
#include "kseq.h"

//KSEQ_INIT(gzFile, gzread)

bwt_seed_t **bwt_init_seed(int seed_n, int max_hit)
{
    int i;
    bwt_seed_t **v = (bwt_seed_t**)malloc(sizeof(bwt_seed_t*));
    (*v) = (bwt_seed_t*)malloc(seed_n * sizeof(bwt_seed_t));
    for (i = 0; i < seed_n; ++i) {
        (*v)[i].loc = (loc_t*)malloc(max_hit * sizeof(loc_t));
        (*v)[i].m = max_hit; (*v)[i].n = 0;
    }
    return v;
}

void bwt_free_seed(bwt_seed_t **v, int seed_n)
{
    int i; 
    for (i = 0; i < seed_n; ++i) free((*v)[i].loc);
    free(*v); free(v);
}

void bwt_set_seed(bwt_seed_t *v, int ref_id, uint8_t is_rev, bwtint_t ref_pos, int read_i)
{
    if (v[read_i].n == v[read_i].m) {
        v[read_i].n = 0;
        return;
    }
    v[read_i].loc[v[read_i].n].ref_id = ref_id;
    v[read_i].loc[v[read_i].n].is_rev = is_rev;
    v[read_i].loc[v[read_i].n].ref_pos = ref_pos;
    ++v[read_i].n;
}

void bwt_init_dp(bwt_seed_t **seed_v, int seed_n)
{
    int i, j;
    for (i = 0; i < seed_n; ++i) {
        for (j = 0; j < (*seed_v)[i].n; ++j) {
            (*seed_v)[i].loc[j].from = START_NODE;
            (*seed_v)[i].loc[j].score = 1;
            (*seed_v)[i].loc[j].NM = 0;
            (*seed_v)[i].loc[j].node_n = 1;
            (*seed_v)[i].loc[j].track_flag = 0;
        }
    }
}

int bwt_seed_con(bwt_seed_t seed1, int l1, bwt_seed_t seed2, int l2, int *con_flag)
{
    if (seed1.loc[l1].is_rev != seed2.loc[l2].is_rev || seed1.loc[l1].ref_id != seed2.loc[l2].ref_id) {
        *con_flag = F_UNCONNECT;
        return 0;
    }
    int dis = (seed1.loc[l1].is_rev?-1:1)*(seed1.loc[l1].ref_pos - seed2.loc[l2].ref_pos) - (seed1.pos - seed2.pos);
    if (abs(dis) <= 1) *con_flag = F_MATCH;
    //if (abs(dis) <= 10) *con_flag = F_MATCH;
    else *con_flag = F_UNCONNECT;
    return abs(dis);
}

void bwt_update_dp(bwt_seed_t **seed_v, int seed_n)
{
    int i, j, m, n, max_dis, con_flag;
    for (i = 1; i < seed_n; ++i) {
        for (j = 0; j <(*seed_v)[i].n; ++j) {
            for (m = i-1; m >= 0; --m) {
                for (n = 0; n < (*seed_v)[m].n; ++n) {
                    max_dis = bwt_seed_con((*seed_v)[i], j, (*seed_v)[m], n, &con_flag);
                    if (con_flag == F_MATCH) {
                        (*seed_v)[i].loc[j].from = (line_node){m, n};
                        (*seed_v)[i].loc[j].score = (*seed_v)[m].loc[n].score+1;
                        (*seed_v)[i].loc[j].NM = (*seed_v)[m].loc[n].NM+max_dis;
                        (*seed_v)[i].loc[j].node_n = (*seed_v)[m].loc[n].node_n +1;
goto Next;
					}
                }
            }
Next:;
        }
    }
}

// return: res_mul_max
int bwt_backtrack(bwt_seed_t **seed_v, reg_t reg, lamsa_aln_para AP, int seed_n, line_node **line, int *node_n)
{
    int l_n = 0; 
    int i, j, k, reg_hit;
    int max_score=0; int track_flag;

    extern node_score *node_init_score(int n);
    extern int heap_add_node(node_score *ns, line_node node, int score, int NM);
    extern void node_free_score(node_score *ns);

    node_score *ns = node_init_score(AP.res_mul_max);

    ns->node_n = 0;
    for (i = seed_n-1; i >= 0; --i) {
        for (j = 0; j < (*seed_v)[i].n; ++j) {
            if ((*seed_v)[i].loc[j].track_flag) continue;
            // check for reg
            reg_hit=0;
            for (k = 0; k < reg.beg_n; ++k) {
                if ((*seed_v)[i].loc[j].ref_id == reg.ref_beg[k].chr-1 && (*seed_v)[i].loc[j].is_rev == reg.ref_beg[k].is_rev
                        && abs((*seed_v)[i].loc[j].ref_pos-reg.ref_beg[k].ref_pos) < AP.SV_len_thd) {
                    reg_hit = 1;
                    break;
                }
            }
            if (reg_hit == 0) {
                for (k = 0; k < reg.end_n; ++k) {
                    if ((*seed_v)[i].loc[j].ref_id == reg.ref_end[k].chr-1 && (*seed_v)[i].loc[j].is_rev == reg.ref_end[k].is_rev
                            && abs((*seed_v)[i].loc[j].ref_pos-reg.ref_end[k].ref_pos) < AP.SV_len_thd) {
                        reg_hit = 1;
                        break;
                    }
                }
            }
            if (reg_hit) (*seed_v)[i].loc[j].score += (*seed_v)[i].loc[j].score/2;
            if ((*seed_v)[i].loc[j].score > max_score) {
                max_score = (*seed_v)[i].loc[j].score;
                heap_add_node(ns, (line_node){i,j}, max_score, (*seed_v)[i].loc[j].NM);
            } else if ((*seed_v)[i].loc[j].score >= max_score/2) 
                heap_add_node(ns, (line_node){i,j}, (*seed_v)[i].loc[j].score, (*seed_v)[i].loc[j].NM);
            // set track_flag
            line_node from = (line_node){i,j};
            while(from.x != START_NODE.x) {
                (*seed_v)[from.x].loc[from.y].track_flag = 1;
                from = (*seed_v)[from.x].loc[from.y].from;
            }
        }
    }
    line_node right; int tmp1, tmp2;
    // track_flag: 0<-untracked(useless), 1<-tracked(candidate), 2<-output(filtered_out)
    while (1) {
        right = node_pop(ns, &tmp1, &tmp2);
        if (right.x == START_NODE.x) break;
        if ((*seed_v)[right.x].loc[right.y].score < max_score/2) continue;
        track_flag = 1;
        line_node from = right;
        while (from.x != START_NODE.x) {
            if ((*seed_v)[from.x].loc[from.y].track_flag == 2) {
                track_flag = 2;
                break;
            }
            from = (*seed_v)[from.x].loc[from.y].from;
        }
        if (track_flag == 2) continue;

        i = (*seed_v)[right.x].loc[right.y].node_n - 1;
        node_n[l_n] = i+1;
        while (right.x != START_NODE.x) {
            line[l_n][i] = right;
            (*seed_v)[right.x].loc[right.y].track_flag = 2;
            right = (*seed_v)[right.x].loc[right.y].from;
            --i;
        }
        l_n++; 
    }

    node_free_score(ns);
    return l_n;
}

int bwt_cluster_seed(bwt_seed_t **seed_v, reg_t reg, lamsa_aln_para AP, int seed_n, line_node **line, int *node_n)
{
    bwt_init_dp(seed_v, seed_n); bwt_update_dp(seed_v, seed_n);
    return bwt_backtrack(seed_v, reg, AP, seed_n, line, node_n);
}

void bwt_set_bound(bwt_seed_t *seed_v, line_node *line, int node_n, int seed_len, int reg_len, bwt_bound *left, bwt_bound *right)
{
    if (seed_v[line[0].x].loc[line[0].y].is_rev) { //rev
        (*right).read_pos = reg_len+1-seed_v[line[0].x].pos;
        (*left).read_pos = reg_len+1-(seed_v[line[node_n-1].x].pos+seed_len-1);

        (*right).ref_pos = seed_v[line[0].x].loc[line[0].y].ref_pos+seed_len-1;
        (*left).ref_pos = seed_v[line[node_n-1].x].loc[line[node_n-1].y].ref_pos;
    } else {
        (*left).read_pos = seed_v[line[0].x].pos;
        (*right).read_pos = seed_v[line[node_n-1].x].pos+seed_len-1;

        (*left).ref_pos = seed_v[line[0].x].loc[line[0].y].ref_pos;
        (*right).ref_pos = seed_v[line[node_n-1].x].loc[line[node_n-1].y].ref_pos+seed_len-1;
    }
}

void bwt_aln_res(int ref_id, uint8_t is_rev, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, uint8_t **read_rbseq, int reg_beg, int reg_len,
                 bwt_bound *left, bwt_bound *right, lamsa_aln_para AP, kseq_t *seqs, line_aln_res *la)
{
    int read_len = seqs->seq.l;
    int extra_ext_len = 100;
    la->cur_res_n = 0;
	la->res[0].chr = ref_id+1; la->res[0].nsrand = 1-is_rev;
    int i; 
    int left_eta_len, extra_beg, right_eta_len, extra_end, extra_reg_len;
    uint8_t *query;
    if (is_rev) {
        left_eta_len = ((read_len-(reg_beg+reg_len-1))>extra_ext_len) ? extra_ext_len:(read_len-(reg_beg+reg_len-1));
        extra_end = reg_beg+reg_len+left_eta_len-1;
        right_eta_len = (reg_beg-1)>extra_ext_len?extra_ext_len:(reg_beg-1);
        extra_beg = reg_beg-right_eta_len;

        extra_reg_len = reg_len + left_eta_len + right_eta_len;
        query = (uint8_t*)malloc(extra_reg_len * sizeof(uint8_t)); // reg_len + extra_len

        for(i = 0; i < extra_reg_len; ++i) 
            query[extra_reg_len-1-i] = (read_bseq[extra_beg-1+i]<4)?3-read_bseq[extra_beg-1+i]:4;
    } else {
        left_eta_len = (reg_beg-1)>extra_ext_len?extra_ext_len:(reg_beg-1);
        extra_beg = reg_beg-left_eta_len;
        right_eta_len = ((read_len-(reg_beg+reg_len-1))>extra_ext_len) ? extra_ext_len:(read_len-(reg_beg+reg_len-1));
        extra_end = reg_beg + reg_len + right_eta_len -1;

        extra_reg_len = reg_len + left_eta_len + right_eta_len;
        query = (uint8_t*)malloc(extra_reg_len * sizeof(uint8_t)); // reg_len + extra_len

        for (i = 0; i < extra_reg_len; ++i) 
            query[i] = read_bseq[extra_beg-1+i];
    }

    uint64_t ref_start = (left->ref_pos - (left->read_pos-1+left_eta_len) - AP.bwt_seed_len < 1) ? 1 : (left->ref_pos - (left->read_pos-1+left_eta_len) - AP.bwt_seed_len);
    int ref_len = (int)(right->ref_pos + reg_len-right->read_pos+right_eta_len + AP.bwt_seed_len - ref_start + 1);
    uint8_t *target = (uint8_t*)malloc(ref_len * sizeof(uint8_t));
    pac2fa_core(bns, pac, ref_id+1, ref_start-1, &ref_len, target);

    int _qle, _tle, qlen, tlen; uint8_t *_q=0, *_t=0;
    qlen = right->read_pos-left->read_pos+1, tlen = right->ref_pos-left->ref_pos+1;
    cigar32_t *mid_cigar; int mid_cigar_n;
    res_t *cur_res = la->res+la->cur_res_n;

    // push head 'S' (before left_extra)
	if (is_rev) _push_cigar1(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), (read_len-extra_end)<<4 | CSOFT_CLIP);
	else _push_cigar1(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), ((extra_beg-1) << 4) | CSOFT_CLIP);

    // mid global-sw
    ksw_global(qlen, query+(left->read_pos-1+left_eta_len), tlen, target+left->ref_pos-ref_start, 5, bwasw_sc_mat, AP.gapo, AP.gape, abs(qlen-tlen)+3, &mid_cigar_n, &mid_cigar);

    cigar32_t *cigar=NULL; int cigar_n, cigar_m;
    if (left->read_pos > 1 || left_eta_len > 0) { // left extend
        // invert query and target,
        qlen = left->read_pos - 1 + left_eta_len;
        _q = (uint8_t*)malloc(qlen * sizeof(uint8_t));
        for (i = 0; i < qlen; ++i) _q[i] = query[qlen-1-i];
        tlen = left->ref_pos-ref_start;
        _t = (uint8_t*)malloc(tlen * sizeof(uint8_t));
        for (i = 0; i < tlen; ++i) _t[i] = target[tlen-1-i];
        // ksw_extend
        ksw_extend_core(qlen, _q, tlen, _t, 5, bwasw_sc_mat, abs(qlen-tlen)+3, AP.bwt_seed_len*bwasw_sc_mat[0], AP, &_qle, &_tle, &cigar, &cigar_n, &cigar_m);
        if (cigar!=NULL) {
            left->read_pos -= _qle;
            left->ref_pos -= _tle;
            _push_cigar1(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), ((qlen-_qle)<<4)|CSOFT_CLIP);
            _invert_cigar(&cigar, cigar_n);
            _push_cigar(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), cigar, cigar_n);
            free(cigar);
        } else _push_cigar1(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), (qlen<<4)|CSOFT_CLIP);
		
        free(_q); free(_t);
    } 
    cur_res->offset = left->ref_pos;
    // push mid-cigar
    _push_cigar(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), mid_cigar, mid_cigar_n);
    free(mid_cigar);

    if (right->read_pos < reg_len || right_eta_len > 0) { // right extend
        qlen = reg_len - right->read_pos + right_eta_len;
        tlen = ref_start+ref_len-1-right->ref_pos;
        ksw_extend_core(qlen, query+right->read_pos+left_eta_len, tlen, target+ref_len-tlen, 5, bwasw_sc_mat, abs(qlen-tlen)+3, AP.bwt_seed_len*bwasw_sc_mat[0], AP, &_qle, &_tle, &cigar, &cigar_n, &cigar_m);
        if (cigar!=NULL) {
            _push_cigar(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), cigar, cigar_n);
            right->read_pos += _qle;
            right->ref_pos += _tle;
			free(cigar);
		}
		_push_cigar1(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), ((reg_len-right->read_pos+right_eta_len)<<4)|CSOFT_CLIP);
    }

    // push tail 'S' (after right_extra)
	if (is_rev)_push_cigar1(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), ((extra_beg-1) << 4) | CSOFT_CLIP);
	else _push_cigar1(&(cur_res->cigar), &(cur_res->cigar_len), &(cur_res->c_m), (read_len-extra_end)<<4 | CSOFT_CLIP);

	if (is_rev) {
        if (*read_rbseq == NULL) {
            *read_rbseq = (uint8_t*)calloc(read_len, sizeof(uint8_t));
            for (i = 0; i < read_len; ++i) (*read_rbseq)[i] = (read_bseq[read_len-1-i]<4)?3-read_bseq[read_len-1-i]:4;
        }
		lamsa_res_aux(la, bns, pac, *read_rbseq, read_len, AP, seqs);
	} else lamsa_res_aux(la, bns, pac, read_bseq, read_len, AP, seqs);

    free(query); free(target);
}

int bwt_aln_core(bwt_t *bwt, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, uint8_t **read_rbseq, reg_t reg, lamsa_aln_para AP, kseq_t *seqs, aln_res *re_res)
{
    int i, j, seed_len = AP.bwt_seed_len, is_rev, ref_id;
    uint8_t *bwt_seed = (uint8_t*)malloc(seed_len * sizeof(uint8_t));
	int max_hit = 100;
    int reg_beg = reg.beg, reg_len = reg.end-reg_beg+1;

    bwt_seed_t **seed_v = bwt_init_seed(reg_len-seed_len+1, max_hit);

    // locate bwt-seeds (exact-match seeds)
    for (i = 0; i <= reg_len-seed_len; ++i) {
        (*seed_v)[i].pos = i+1; // 1-base, on reg
        for (j = i; j < i+seed_len; ++j) bwt_seed[j-i] = read_bseq[reg_beg-1+j];
        uint64_t k = 0, l = bwt->seq_len, m;
		if (!bwt_match_exact_alt(bwt, seed_len, bwt_seed, &k, &l)) continue;
        if (l-k+1 <= max_hit) { // set_seed directly
            for (m = k; m <=l; ++m) {
                bwtint_t ref_pos = bwt_sa(bwt, m);
                bwtint_t pos = bns_depos(bns, ref_pos, &is_rev);
                bns_cnt_ambi(bns, pos, seed_len, &ref_id);
                bwt_set_seed(*seed_v, ref_id, is_rev, pos-bns->anns[ref_id].offset+1-(is_rev?(seed_len-1):0), i);
            }
        } else if (l-k+1 <= 5 * max_hit){ 
			// select specific seed-results, which are close to the existing result.
			// keep 100 seed-results, at most.
            int cnt=0, reg_hit;
			for (m = k; m <=l; ++m) {
                bwtint_t ref_pos = bwt_sa(bwt, m);
                bwtint_t pos = bns_depos(bns, ref_pos, &is_rev);
                bns_cnt_ambi(bns, pos, seed_len, &ref_id);
                uint64_t abs_pos = pos-bns->anns[ref_id].offset+1-(is_rev?(seed_len-1):0);
                reg_hit = 0;
                for (j = 0; j < reg.beg_n; ++j) {
                    if (ref_id == reg.ref_beg[j].chr-1 && is_rev == reg.ref_beg[j].is_rev
                     && abs(abs_pos-reg.ref_beg[j].ref_pos) < AP.SV_len_thd) {
                        reg_hit = 1;
                        break;
                    }
                }
                if (reg_hit == 0) {
                    for (j = 0; j < reg.end_n; ++j) {
                        if (ref_id == reg.ref_end[j].chr-1 && is_rev == reg.ref_end[j].is_rev
                                && abs(abs_pos-reg.ref_end[j].ref_pos) < AP.SV_len_thd) {
                            reg_hit = 1;
                            break;
                        }
                    }
                }
                if (reg_hit == 0) continue;
                bwt_set_seed(*seed_v, ref_id, is_rev, abs_pos, i);
                cnt++;
                if (cnt == max_hit) break;
            }
        }
    }
    // cluster seeds, find optimal res_mul_max cluster
    line_node **line = (line_node**)malloc(AP.res_mul_max * sizeof(line_node*));
    for (i = 0; i < AP.res_mul_max; ++i) line[i] = (line_node*)malloc((reg_len-seed_len+1) * sizeof(line_node));
    int *node_n = (int*)malloc(AP.res_mul_max * sizeof(int)); int l_n;

    bwt_bound left_bound, right_bound;
    // use DP, Spanning-tree and Pruning
    extern void aln_reloc_res(aln_res *a_res, int line_n, int XA_m);
    if ((l_n = bwt_cluster_seed(seed_v, reg, AP, reg_len-seed_len+1, line, node_n)) > 0) {
        for (i = 0; i < l_n; ++i) {
            if (re_res->l_n == re_res->l_m) aln_reloc_res(re_res, (re_res->l_n << 1), AP.res_mul_max);

            bwt_set_bound(*seed_v, line[i], node_n[i], seed_len, reg_len, &left_bound, &right_bound);
            bwt_aln_res((*seed_v)[line[i][0].x].loc[line[i][0].y].ref_id, (*seed_v)[line[i][0].x].loc[line[i][0].y].is_rev, bns, pac, read_bseq, read_rbseq, reg_beg, reg_len, &left_bound, &right_bound, AP, seqs, re_res->la+re_res->l_n);
            re_res->la[re_res->l_n].line_score = 0;
            re_res->l_n++;
        }
    }
    // free
    free(bwt_seed); bwt_free_seed(seed_v, reg_len-seed_len+1);
    for (i = 0; i < AP.res_mul_max; ++i) free(line[i]); free(line); free(node_n);
    return 0;
}

void bwt_aln_remain(aln_reg *a_reg, aln_res *re_res, bwt_t *bwt, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, uint8_t **read_rbseq, lamsa_aln_para AP, kseq_t *seqs)
{
    aln_reg *re_reg = aln_init_reg(seqs->seq.l); 
    if (get_remain_reg(a_reg, re_reg, AP, AP.bwt_max_len) == 0) goto End;

    int i;
    re_res->l_n = 0;
    for (i = 0; i < re_reg->reg_n; ++i)
        bwt_aln_core(bwt, bns, pac, read_bseq, read_rbseq, re_reg->reg[i], AP, seqs, re_res);

End:
    aln_free_reg(re_reg);
}
