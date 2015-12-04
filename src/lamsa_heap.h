#ifndef LAMSA_HEAP
#define LAMSA_HEAP


line_node node_pop(node_score *ns, int *score, int *NM);
//void node_max_heap(node_score *ns, int i);
void build_node_max_heap(node_score *ns);
line_node node_heap_extract_max(node_score *ns, int *score);
//void build_node_maxpos_heap(node_score *ns);
//line_node node_heap_extract_maxpos(node_score *ns);
void build_node_minpos_heap(node_score *ns);
line_node node_heap_extract_minpos(node_score *ns);

//void node_min_heap(node_score *ns, int i);
void build_node_min_heap(node_score *ns);
//line_node node_heap_extract_min(node_score *ns);
int node_heap_update_min(node_score *ns, line_node node, int score, int NM);


#endif
