#ifndef LAMSA_DP_CON_H
#define LAMSA_DP_CON_H

#include "lamsa_aln.h"
#include "frag_check.h"

int frag_line_remain(aln_reg *a_reg, aln_msg *a_msg, frag_msg **f_msg, lamsa_aln_per_para APP, lamsa_aln_para AP, kseq_t *seq,
                     line_node **line, int *line_end, int *line_m, frag_dp_node ***f_node, line_node **_line, int *_line_end,
                     int line_n_max);
int frag_line_BCC(aln_msg *a_msg, frag_msg **f_msg,
                  lamsa_aln_per_para APP, lamsa_aln_para AP, kseq_t *seq,
                  line_node **line, int *line_end, int *line_m,
                  frag_dp_node ***f_node,
                  line_node **_line, int line_n_max);

#endif
