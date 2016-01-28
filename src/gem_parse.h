#ifndef GEM_PARSE_H
#define GEM_PARSE_H
#include "lamsa_aln.h"

int gem_map_read(FILE *mapf, map_msg *m_msg, char *gem_line, int line_size);
int gem_map_msg(map_msg *m_msg, int max_n);
map_msg *map_init_msg(int n);
void map_free_msg(map_msg *m_msg, int n);

#endif
