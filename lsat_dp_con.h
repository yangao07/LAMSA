#ifndef LSAT_DP_CON_H
#define LSAT_DP_CON_H

#include "lsat_aln.h"
#include "frag_check.h"

#define S_REMAIN_CON 3

int frag_line_remain(aln_reg *a_reg, aln_msg *a_msg, frag_msg **f_msg, lsat_aln_per_para APP, lsat_aln_para AP, 
                     line_node **line, int *line_end, int *line_m, frag_dp_node ***f_node, line_node **_line, int *_line_end,
                     int line_n_max, int per_max_multi);
int frag_line_BCC(aln_msg *a_msg, frag_msg **f_msg,
                  lsat_aln_per_para APP, lsat_aln_para AP,
                  line_node **line, int *line_end, int *line_m,
                  frag_dp_node ***f_node,
                  line_node **_line, int line_n_max, int per_max_multi);

#endif
