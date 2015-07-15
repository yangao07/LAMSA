#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "frag_check.h"
#include "lsat_aln.h"
#include "split_mapping.h"
#include "bwt_aln.h"
#include "bntseq.h"
#include "bwt.h"
#include "ksw.h"

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
    if (abs(dis) <= 10) *con_flag = F_MATCH;
    else *con_flag = F_UNCONNECT;
    return abs(dis);
}

void bwt_update_dp(bwt_seed_t **seed_v, int seed_n)
{
    int i, j, m, n, max_score=0, max_dis, con_flag;
    line_node max_from;
    for (i = 1; i < seed_n; ++i) {
        for (j = 0; j <(*seed_v)[i].n; ++j) {
            for (m = i-1; m >= 0; --m) {
                for (n = 0; n < (*seed_v)[m].n; ++n) {
                    max_dis = bwt_seed_con((*seed_v)[i], j, (*seed_v)[m], n, &con_flag);
                    if (con_flag == F_MATCH) {
                        (*seed_v)[i].loc[j].from = (line_node){m, n};
                        (*seed_v)[i].loc[j].score = (*seed_v)[m].loc[n].score+1;
goto Next;
					}
                }
            }
Next:;
        }
    }
}

int bwt_backtrack(bwt_seed_t **seed_v, int seed_n, line_node **line)
{
    int i, j;
    int max_score=0/*also is max_node_num*/; line_node max_node=START_NODE;
    for (i = seed_n-1; i >= 0; --i) {
        for (j = 0; j < (*seed_v)[i].n; ++j) {
            if ((*seed_v)[i].loc[j].score > max_score) {
                max_score = (*seed_v)[i].loc[j].score;
                max_node = (line_node){i, j};
            }
        }
    }
    line_node right = max_node; i = max_score - 1;
    while (right.x != START_NODE.x) {
        (*line)[i] = right; 
        right = (*seed_v)[right.x].loc[right.y].from;
        --i;
    }
    return max_score;
}

int bwt_cluster_seed(bwt_seed_t **seed_v, int seed_n, line_node **line)
{
    bwt_init_dp(seed_v, seed_n); bwt_update_dp(seed_v, seed_n);
    return bwt_backtrack(seed_v, seed_n, line);
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
                 bwt_bound *left, bwt_bound *right, lsat_aln_para AP, lsat_aln_per_para APP, line_aln_res *la)
{
    la->cur_res_n = 0;
	la->res[0].chr = ref_id+1; la->res[0].nsrand = 1-is_rev;
    int i; uint8_t *query = (uint8_t*)malloc(reg_len * sizeof(uint8_t));
    if (is_rev) for (i = 0; i < reg_len; ++i) query[reg_len-1-i] = (read_bseq[reg_beg-1+i]<4)?3-read_bseq[reg_beg-1+i]:4; 
    else for (i = 0; i < reg_len; ++i) query[i] = read_bseq[reg_beg-1+i];

    uint64_t ref_start = (left->ref_pos - (left->read_pos - 1) - AP.bwt_seed_len < 1) ? 1 : (left->ref_pos - (left->read_pos - 1) - AP.bwt_seed_len);
    int ref_len = (int)(right->ref_pos + reg_len - right->read_pos + AP.bwt_seed_len - ref_start + 1);
    uint8_t *target = (uint8_t*)malloc(ref_len * sizeof(uint8_t));
    pac2fa_core(bns, pac, ref_id+1, ref_start-1, &ref_len, target);

    int _qle, _tle, qlen, tlen; uint8_t *_q=0, *_t=0;
    qlen = right->read_pos-left->read_pos+1, tlen = right->ref_pos-left->ref_pos+1;
    cigar32_t *mid_cigar; int mid_cigar_n;
    // push head 'S'
	if (is_rev) _push_cigar1(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), (APP.read_len-(reg_beg+reg_len-1))<<4 | CSOFT_CLIP);
	else _push_cigar1(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), ((reg_beg-1) << 4) | CSOFT_CLIP);

    // mid global-sw
    ksw_global(qlen, query+left->read_pos-1, tlen, target+left->ref_pos-ref_start, 5, bwasw_sc_mat, AP.gapo, AP.gape, abs(qlen-tlen)+3, &mid_cigar_n, &mid_cigar);

    cigar32_t *cigar=NULL; int cigar_n, cigar_m;
    if (left->read_pos > 1) { // left extend
        // invert query and target,
        qlen = left->read_pos-1;
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
            _push_cigar1(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), ((left->read_pos-1)<<4)|CSOFT_CLIP);
            _invert_cigar(&cigar, cigar_n);
            _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), cigar, cigar_n);
            free(cigar);
        } else _push_cigar1(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), ((left->read_pos-1)<<4)|CSOFT_CLIP);
		
        free(_q); free(_t);
    } 
    la->res[la->cur_res_n].offset = left->ref_pos;
    // push mid-cigar
    _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), mid_cigar, mid_cigar_n);
    free(mid_cigar);

    if (right->read_pos < reg_len) { // right extend
        qlen = reg_len - right->read_pos;
        tlen = ref_start+ref_len-1-right->ref_pos;
        ksw_extend_core(qlen, query+right->read_pos, tlen, target+ref_len-tlen, 5, bwasw_sc_mat, abs(qlen-tlen)+3, AP.bwt_seed_len*bwasw_sc_mat[0], AP, &_qle, &_tle, &cigar, &cigar_n, &cigar_m);
        if (cigar!=NULL) {
            _push_cigar(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), cigar, cigar_n);
            right->read_pos += _qle;
            right->ref_pos += _tle;
			free(cigar);
		}
		_push_cigar1(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), ((reg_len-right->read_pos)<<4)|CSOFT_CLIP);
    }

    // push tail 'S'
	if (is_rev)_push_cigar1(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), ((reg_beg-1) << 4) | CSOFT_CLIP);
	else _push_cigar1(&(la->res[la->cur_res_n].cigar), &(la->res[la->cur_res_n].cigar_len), &(la->res[la->cur_res_n].c_m), (APP.read_len-(reg_beg+reg_len-1))<<4 | CSOFT_CLIP);

	if (is_rev) {
        if (*read_rbseq == NULL) {
            *read_rbseq = (uint8_t*)calloc(APP.read_len, sizeof(uint8_t));
            for (i = 0; i < APP.read_len; ++i) (*read_rbseq)[i] = (read_bseq[APP.read_len-1-i]<4)?3-read_bseq[APP.read_len-1-i]:4;
        }
		lsat_res_aux(la, bns, pac, *read_rbseq, APP.read_len, AP, APP);
	}
	else lsat_res_aux(la, bns, pac, read_bseq, APP.read_len, AP, APP);

    free(query); free(target);
}

int bwt_aln_core(bwt_t *bwt, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, uint8_t **read_rbseq, reg_t reg, lsat_aln_para AP, lsat_aln_per_para APP, aln_res *re_res)
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
            int cnt=0;
			for (m = k; m <=l; ++m) {
                bwtint_t ref_pos = bwt_sa(bwt, m);
                bwtint_t pos = bns_depos(bns, ref_pos, &is_rev);
				if (is_rev != reg.is_rev) continue;
                bns_cnt_ambi(bns, pos, seed_len, &ref_id);
				if (ref_id != reg.refid) continue;
                uint64_t abs_pos = pos-bns->anns[ref_id].offset+1-(is_rev?(seed_len-1):0);
                //bwt_set_seed(*seed_v, ref_id, is_rev, abs_pos, i);
				if (is_rev) {
					//fprintf(stderr, "ref_beg: %d\nabs_pos: %d\nref_end: %d\n", reg.ref_beg, abs_pos, reg.ref_end);
					if ((!reg.ref_beg || (abs_pos < reg.ref_beg && reg.ref_beg-abs_pos < AP.SV_len_thd)) && (!reg.ref_end || (abs_pos > reg.ref_end && abs_pos-reg.ref_end < AP.SV_len_thd))) {
						bwt_set_seed(*seed_v, ref_id, is_rev, abs_pos, i);
						cnt++;
						if (cnt == max_hit) break;
					}
				} else {
					if ((reg.ref_beg==0 || (abs_pos > reg.ref_beg && abs_pos-reg.ref_beg<AP.SV_len_thd)) && (!reg.ref_end || (abs_pos < reg.ref_end && reg.ref_end-abs_pos < AP.SV_len_thd))) {
						bwt_set_seed(*seed_v, ref_id, is_rev, abs_pos, i);
						cnt++;
						if (cnt == max_hit) break;
					}
				}
			}
		}
    }
    // cluster seeds, find one best cluster //XXX
    line_node *line = (line_node*)malloc((reg_len-seed_len+1) * sizeof(line_node));
    int node_n; bwt_bound left_bound, right_bound;
	// XXX use DP, Spanning-tree and Pruning
    if ((node_n = bwt_cluster_seed(seed_v, reg_len-seed_len+1, &line)) > 0) {
		if (re_res->l_n == re_res->l_m) {
			re_res->l_m <<= 1;
			if ((re_res->la = (line_aln_res*)realloc(re_res->la, re_res->l_m * sizeof(line_aln_res))) == NULL) { fprintf(stderr, "[bwt_aln_core] Not enough memory for aln_res.\n"); exit(0); }
			for (i = (re_res->l_m)>>1; i < re_res->l_m; ++i) {
				re_res->la[i].res_m = 10, re_res->la[i].cur_res_n = 0, re_res->la[i].split_flag = 0;
				re_res->la[i].res = (res_t*)malloc(10 * sizeof(res_t));
				for (j = 0; j < 10; ++j) {
					re_res->la[i].res[j].c_m = 100;
					re_res->la[i].res[j].cigar = (cigar32_t*)malloc(100 * sizeof(cigar32_t));
					re_res->la[i].res[j].cigar_len = 0;
				}
				re_res->la[i].tol_score = re_res->la[i].tol_NM = 0;
				re_res->la[i].trg_m = 10, re_res->la[i].trg_n = 0;
				re_res->la[i].trg = (line_node*)malloc(10 * sizeof(line_node));
			}
		}
		bwt_set_bound(*seed_v, line, node_n, seed_len, reg_len, &left_bound, &right_bound);
		bwt_aln_res((*seed_v)[line[0].x].loc[line[0].y].ref_id, (*seed_v)[line[0].x].loc[line[0].y].is_rev, bns, pac, read_bseq, read_rbseq, reg_beg, reg_len, &left_bound, &right_bound, AP, APP, re_res->la+re_res->l_n);
		re_res->la[re_res->l_n].merg_msg = (line_node){1,-1};
        re_res->l_n++;
    }
    // free
    free(bwt_seed); bwt_free_seed(seed_v, reg_len-seed_len+1);
    free(line);
    return 0;
}

void bwt_aln_remain(aln_reg *a_reg, aln_res *re_res, bwt_t *bwt, bntseq_t *bns, uint8_t *pac, uint8_t *read_bseq, uint8_t **read_rbseq, lsat_aln_per_para APP, lsat_aln_para AP)
{
    aln_reg *re_reg = aln_init_reg(APP.read_len); 
    if (get_remain_reg(a_reg, re_reg, AP, 300) == 0) goto End; //XXX

    int i;
    re_res->l_n = 0;
    // extend the remain_reg or not?XXX
    for (i = 0; i < re_reg->reg_n; ++i)
        bwt_aln_core(bwt, bns, pac, read_bseq, read_rbseq, re_reg->reg[i], AP, APP, re_res);

End:
    aln_free_reg(re_reg);
}
