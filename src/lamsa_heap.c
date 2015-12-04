#include <stdlib.h>
#include <stdio.h>
#include "lamsa_aln.h"

line_node node_pop(node_score *ns, int *score, int *NM)
{
    if (ns->node_n < 1) 
        return (line_node){-1,0};
    line_node node = ns->node[--(ns->node_n)];
    *score = ns->score[ns->node_n];
    *NM = ns->NM[ns->node_n];
    return node;
}

//max score
void node_max_heap(node_score *ns, int i)
{
    int l = 2*i+1; int r = 2*(i+1);
    int max, tmp;
    if (l < ns->node_n && (ns->score[l] > ns->score[i] || (ns->score[l] == ns->score[i] && ns->NM[l] < ns->NM[i]))) 
        max = l;
    else max = i;
    if (r < ns->node_n && (ns->score[r] > ns->score[max] || (ns->score[r] == ns->score[max] && ns->NM[r] < ns->NM[max])))
        max = r;
    if (max != i)
    {
        tmp = ns->node[i].x, ns->node[i].x = ns->node[max].x, ns->node[max].x = tmp;
        tmp = ns->node[i].y, ns->node[i].y = ns->node[max].y, ns->node[max].y = tmp;
        tmp = ns->score[i], ns->score[i] = ns->score[max], ns->score[max] = tmp;
        node_max_heap(ns, max);
    }
}

void build_node_max_heap(node_score *ns)
{
    if (ns->node_n == 0) return;
    int i;
    for (i = (ns->node_n-1)/2; i >= 0; --i)
        node_max_heap(ns, i);
}

line_node node_heap_extract_max(node_score *ns, int *score)
{
    if (ns->node_n < 1) {
        //fprintf(stderr, "[heap_sort] MAX error.\n"); exit(1);
        return (line_node){-1,0};
    }
    line_node max = ns->node[0];
    *score = ns->score[0];
    ns->node[0] = ns->node[--(ns->node_n)];
    ns->score[0] = ns->score[ns->node_n];
    ns->NM[0] = ns->NM[ns->node_n];
    node_max_heap(ns, 0);
    return max;
}

// max end-pos
void node_maxpos_heap(node_score *ns, int i)
{
    int l = 2*i+1; int r = 2*(i+1);
    int max, tmp;
    if (l < ns->node_n && ns->node[l].x > ns->node[i].x) 
        max = l;
    else max = i;
    if (r < ns->node_n && ns->node[r].x > ns->node[max].x) 
        max = r;
    if (max != i)
    {
        tmp = ns->node[i].x, ns->node[i].x = ns->node[max].x, ns->node[max].x = tmp;
        tmp = ns->node[i].y, ns->node[i].y = ns->node[max].y, ns->node[max].y = tmp;
        tmp = ns->score[i], ns->score[i] = ns->score[max], ns->score[max] = tmp;
        tmp = ns->NM[i], ns->NM[i] = ns->NM[max], ns->NM[max] = tmp;
        node_maxpos_heap(ns, max);
    }
}

void build_node_maxpos_heap(node_score *ns)
{
    int i;
    for (i = (ns->node_n-1)/2; i >= 0; --i)
        node_maxpos_heap(ns, i);
}

line_node node_heap_extract_maxpos(node_score *ns)
{
    if (ns->node_n < 1)
    {
        //fprintf(stderr, "[heap_sort] MAX error.\n"); exit(1);
        return (line_node){-1,0};
    }
    line_node max = ns->node[0];
    ns->node[0] = ns->node[--(ns->node_n)];
    ns->score[0] = ns->score[ns->node_n];
    ns->NM[0] = ns->NM[ns->node_n];
    node_maxpos_heap(ns, 0);
    return max;
}

// min end-pos
void node_minpos_heap(node_score *ns, int i)
{
    int l = 2*i+1; int r = 2*(i+1);
    int min, tmp;
    if (l < ns->node_n && ns->node[l].x < ns->node[i].x) 
        min = l;
    else min = i;
    if (r < ns->node_n && ns->node[r].x < ns->node[min].x) 
        min = r;
    if (min != i)
    {
        tmp = ns->node[i].x, ns->node[i].x = ns->node[min].x, ns->node[min].x = tmp;
        tmp = ns->node[i].y, ns->node[i].y = ns->node[min].y, ns->node[min].y = tmp;
        tmp = ns->score[i], ns->score[i] = ns->score[min], ns->score[min] = tmp;
        tmp = ns->NM[i], ns->NM[i] = ns->NM[min], ns->NM[min] = tmp;
        node_minpos_heap(ns, min);
    }
}

void build_node_minpos_heap(node_score *ns)
{
    int i;
    for (i = (ns->node_n-1)/2; i >= 0; --i)
        node_minpos_heap(ns, i);
}

line_node node_heap_extract_minpos(node_score *ns)
{
    if (ns->node_n < 1)
    {
        //fprintf(stderr, "[heap_sort] MAX error.\n"); exit(1);
        return (line_node){-1,0};
    }
    line_node min = ns->node[0];
    ns->node[0] = ns->node[--(ns->node_n)];
    ns->score[0] = ns->score[ns->node_n];
    ns->NM[0] = ns->NM[ns->node_n];
    node_minpos_heap(ns, 0);
    return min;
}

/*void node_heap_update_max(node_score *ns, line_node node, int score)
{
    if (ns->score[0] < score)
    {
        ns->score[0] = score;
        ns->node[0] = node;
        node_max_heap(ns, 0);
    }
}*/

void node_min_heap(node_score *ns, int i)
{
    int l = 2*i+1; int r = 2*(i+1);
    int min, tmp;
    if (l < ns->node_n && (ns->score[l] < ns->score[i] || (ns->score[l] == ns->score[i] && ns->NM[l] > ns->NM[i])))
        min = l;
    else min = i;
    if (r < ns->node_n && (ns->score[r] < ns->score[min] || (ns->score[r] == ns->score[min] && ns->NM[r] > ns->NM[min])))
        min = r;
    if (min != i)
    {
        tmp = ns->node[i].x, ns->node[i].x = ns->node[min].x, ns->node[min].x = tmp;
        tmp = ns->node[i].y, ns->node[i].y = ns->node[min].y, ns->node[min].y = tmp;
        tmp = ns->score[i], ns->score[i] = ns->score[min], ns->score[min] = tmp;
        tmp = ns->NM[i], ns->NM[i] = ns->NM[min], ns->NM[min] = tmp;
        node_min_heap(ns, min);
    }
}

void build_node_min_heap(node_score *ns)
{
    int i;
    for (i = (ns->node_n-1)/2; i >= 0; --i)
        node_min_heap(ns, i);
}

/*line_node node_heap_extract_min(node_score *ns)
{
    if (ns->node_n < 1)
    {
        fprintf(stderr, "[heap_sort] MIN error.\n"); exit(1);
    }
    line_node min = ns->node[0];
    ns->node[0] = ns->node[(ns->node_n)--];

    node_min_heap(ns, 0);
    return min;
}*/

int node_heap_update_min(node_score *ns, line_node node, int score, int NM)
{
    if (ns->score[0] < score || (ns->score[0] == score && ns->NM[0] > NM))
    {
        int ret = ns->node[0].x;
        ns->score[0] = score;
        ns->NM[0] = NM;
        ns->node[0] = node;
        node_min_heap(ns, 0);
        return ret;
    } else return -2;
}

#define max_heap(type_t)                \
    void max_heap(type_t *A, int i)     \
    {                                   \
        A[i] = A[0];                    \
    }                                   \

#define HEAP_SORT(type_t) max_heap(type_t)

