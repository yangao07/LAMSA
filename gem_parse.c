#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "gem_parse.h"
#include "lsat_aln.h"
#include "frag_check.h"

#define LINE_SIZE 10000

void md2m(char *c, int *m, int *mm)
{
    *m=0; *mm=0;
    if (*c == 0) return;
    int i, c_i, _c_i;
    long x;
    char *s, *t;

    int len = strlen(c);
    for (i = 0; i < len; ++i) {
        if (c[i] >= 'A' && c[i] <= 'T') //XXX
        {
            (*mm)++;
            c[i] = ' ';
        }
    }
    s = strtok(c, " ");
    (*m) += atoi(s);
    
    while (s != NULL) {
        s = strtok(NULL, " ");
        if (s == NULL) break;
        (*m) += atoi(s);
    }
}

//9A50G11>2-1
void md2cigar(char *md, map_t *map)
{
    int md_len = strlen(md);
    int i; char *c; int c_i, _c_i;
    char indel_c[10], id, m_c[100]; int indel_n;
    c=strtok(md, ">");
    c_i = strlen(c)+1;
    int m, mm;
    md2m(c, &m, &mm);
    map->cigar->cigar_n = 0;
    _push_cigar1(&(map->cigar->cigar), &(map->cigar->cigar_n), &(map->cigar->cigar_m), (m+mm) << 4 | CMATCH);
    map->NM = mm;
    while (c != NULL) {
        if (c_i >= md_len) break;
        c = strtok(md+c_i, ">");
        if (c == NULL) break;
        c_i += (strlen(c)+1);
        _c_i = c_i;
        sscanf(c, "%[^+-]%[+-]%s", indel_c, &id, m_c);
        c_i = _c_i;
        indel_n = atoi(indel_c);
        _push_cigar1(&(map->cigar->cigar), &(map->cigar->cigar_n), &(map->cigar->cigar_m), (indel_n) << 4 | (id=='+'?CDEL:CINS));
        md2m(m_c, &m, &mm);
        _push_cigar1(&(map->cigar->cigar), &(map->cigar->cigar_n), &(map->cigar->cigar_m), (m+mm) << 4 | (CMATCH));
        map->NM += (indel_n+mm);
    }
}

map_msg *map_init_msg(void)
{
    int i, max_n = 1;
    map_msg *m_msg = (map_msg*)malloc(sizeof(map_msg));
    m_msg->map_n = 0; m_msg->map_m = max_n; m_msg->map = (map_t*)malloc(max_n * sizeof(map_t));
    for (i = 0; i < max_n; ++i) {
        m_msg->map[i].cigar = (cigar_t*)malloc(sizeof(cigar_t));
        m_msg->map[i].cigar->cigar_n = 0; m_msg->map[i].cigar->cigar_m = 10; m_msg->map[i].cigar->cigar = (cigar32_t*)malloc(10 * sizeof(cigar32_t));
    }
    return m_msg;
}

void map_free_msg(map_msg *m_msg)
{
    int i;
    for (i = 0; i < m_msg->map_m; ++i) {
        free(m_msg->map[i].cigar->cigar);
        free(m_msg->map[i].cigar);
    }
    free(m_msg->map); free(m_msg);
}

#ifdef _GEM_MAIN_
int main(int argc, char* argv[])
{
    int i; char file[100]; FILE *mapf;
    strcpy(file, argv[1]);
    mapf = fopen(file, "r");
    map_msg *m_msg = (map_msg*)malloc(sizeof(map_msg)); int max_n = 100;
    m_msg->map_n = 0; m_msg->map_m = max_n; m_msg->map = (map_t*)malloc(max_n * sizeof(map_t));
    for (i = 0; i < max_n; ++i) {
        m_msg->map[i].cigar = (cigar_t*)malloc(sizeof(cigar_t));
        m_msg->map[i].cigar->cigar_n = 0; m_msg->map[i].cigar->cigar_m = 10; m_msg->map[i].cigar->cigar = (cigar32_t*)malloc(10 * sizeof(cigar32_t));
    }

    m_msg->map_n = 0;
    char gem_line[10000];
    char name[100], aln_msg[10000];
    char *t, *s; int t_i, _t_i;
    char chr[10], strand, md[200];
    long long offset;


    if (fgets(gem_line, 10000, mapf) == NULL) return -1;
    sscanf(gem_line, "%[^\t]\t%*[^\t]\t%*[^\t]\t%[^\n]\n", name, aln_msg);
    t = strtok(aln_msg, ",");
    t_i = strlen(t)+1;
    while (t != NULL) {
        _t_i = t_i;
        sscanf(t, "%[^:]:%c:%lld:%[^:]", chr, &strand, &offset, md);
        t_i = _t_i;
        md2cigar(md, &(m_msg->map[m_msg->map_n]));
        m_msg->map[m_msg->map_n].offset = offset;
        m_msg->map[m_msg->map_n].strand = strand;
        strcpy(m_msg->map[m_msg->map_n].chr, chr);

        printf("%s, %c, %lld, ", chr, strand, offset);
        printcigar(stdout, m_msg->map[m_msg->map_n].cigar->cigar, m_msg->map[m_msg->map_n].cigar->cigar_n);
        ++(m_msg->map_n);
        t = strtok(aln_msg+t_i, ",");
        if (t == NULL) break;
        t_i += strlen(t)+1;
    }
    map_free_msg(m_msg);
    return 0;
}

#endif

int gem_map_read(FILE *mapf, map_msg *m_msg, int max_n)
{
    m_msg->map_n = 0;
    char gem_line[10000];
    char name[100], aln_msg[10000]; int msg_len;
    char *t, *s; int t_i, _t_i;
    char chr[10], strand, md[200];
    long long offset;
    int i;

    if (fgets(gem_line, 10000, mapf) == NULL) return -1;
    sscanf(gem_line, "%[^\t]\t%*[^\t]\t%*[^\t]\t%[^\n]\n", name, aln_msg);
    if (strcmp(aln_msg, "-") == 0) return 0;
    t = strtok(aln_msg, ",");
    msg_len = strlen(aln_msg);
    t_i = strlen(t)+1;
    while (t != NULL) {
        if (m_msg->map_n == m_msg->map_m) {
            m_msg->map_m <<= 1;
            m_msg->map = (map_t*)realloc(m_msg->map, (m_msg->map_m>>1) * sizeof(map_t));
            if (m_msg->map == NULL) {
                fprintf(stderr, "[gem_map_read] Error: not enough memory.\n");
                exit(0);
            }
            for (i = m_msg->map_m>>1; i < m_msg->map_m; ++i) {
                m_msg->map[i].cigar = (cigar_t*)malloc(sizeof(cigar_t));
                m_msg->map[i].cigar->cigar_n = 0; m_msg->map[i].cigar->cigar_m = 10;
                m_msg->map[i].cigar->cigar = (cigar32_t*)malloc(10 * sizeof(cigar32_t));
            }
        }
        _t_i = t_i;
        sscanf(t, "%[^:]:%c:%lld:%[^:]", chr, &strand, &offset, md);
        t_i = _t_i;
        md2cigar(md, &(m_msg->map[m_msg->map_n]));
        m_msg->map[m_msg->map_n].offset = offset;
        m_msg->map[m_msg->map_n].strand = strand;
        strcpy(m_msg->map[m_msg->map_n].chr, chr);

        //printf("%s, %c, %lld, ", chr, strand, offset);
        //printcigar(stdout, m_msg->map[m_msg->map_n].cigar->cigar, m_msg->map[m_msg->map_n].cigar->cigar_n);
        ++(m_msg->map_n);
        if (t_i >= msg_len) break;
        t = strtok(aln_msg+t_i, ",");
        if (t == NULL) break;
        t_i += strlen(t)+1;
    }

    return m_msg->map_n;
}
