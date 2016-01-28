#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "gem_parse.h"
#ifndef _GEM_MAIN_
#include "lamsa_aln.h"
#include "frag_check.h"
#endif

#ifdef _GEM_MAIN_

void printcigar(FILE *outp, cigar32_t *cigar, int cigar_len)
{
	int i; 
    for (i = 0; i < cigar_len; i++)
		fprintf(outp, "%d%c", (int)(cigar[i]>>4), CIGAR_STR[(int)(cigar[i] & 0xf)]);
	fprintf(stdout, "\n");
}

void _push_cigar1(cigar32_t **cigar, int *cigar_n, int *cigar_m, cigar32_t _cigar) {
	if (_cigar >> 4 == 0) return;
    int i;
    i = *cigar_n;
    if (((i-1) >=0) && (((*cigar)[i-1] & 0xf) == (_cigar & 0xf)))
        (*cigar)[i-1] += ((_cigar >> 4) << 4);
    else
    {
        if (i == *cigar_m) {
            (*cigar_m) <<= 1;
			(*cigar) = (cigar32_t*)realloc(*cigar, (*cigar_m) * sizeof (cigar32_t));
			if ((*cigar) == NULL)	{fprintf(stderr, "\n[frag_check] Memory is not enougy.\n");exit(1);}
        }
        (*cigar)[i] = _cigar;
        ++(*cigar_n);
    }
}
#endif 

void md2m(char *c, int len, int *m, int *mm)
{
    *m=0; *mm=0;
    if (*c == 0) return;
    int i; char *s, *t;

    for (i = 0; i < len; ++i) {
        if (c[i] >= 'A' && c[i] <= 'T')
        {
            (*mm)++;
            c[i] = ' ';
        }
    }
    s = strtok_r(c, " ", &t);
    while (s != NULL) {
        (*m) += atoi(s);
        s = strtok_r(NULL, " ", &t);
    }
}

int get_d(int indel_n)
{
    int i = 1;
    while (indel_n >= 10) {
        i++;
        indel_n /= 10;
    }
    return i;
}

//9A50G11>2-1
void md2cigar(char *md, map_t *map)
{
    int md_len = strlen(md);
    int i=0; 
	char *c = (char*)malloc((md_len+1) * sizeof(char)); int c_i;
	char *indel = (char*)malloc((md_len+1) * sizeof(char));
    char id; int indel_n=0, bd=0, bi=0;
    int m, mm;
	map->cigar->cigar_n = 0;
	map->NM = 0;
	while (i < md_len) {
		if (md[i] == '>') { // indel
			sscanf(md+i, ">%[^+-]%[+-]", indel, &id);
            indel_n = atoi(indel);
            if (id=='+') {
                bd += indel_n;
                _push_cigar1(&(map->cigar->cigar), &(map->cigar->cigar_n), &(map->cigar->cigar_m), (indel_n) << 4 | CDEL);
            } else { 
                bi += indel_n;
                _push_cigar1(&(map->cigar->cigar), &(map->cigar->cigar_n), &(map->cigar->cigar_m), (indel_n) << 4 | CINS);
            }
			map->NM += indel_n;
            i += (get_d(indel_n) + 2);
		} else {
			c_i = i;
			for (; md[i] && md[i]!='>'; ++i) {
				c[i-c_i] = md[i];
			}
            m = 0, mm = 0;
            c[i-c_i] = '\0';
			md2m(c, i-c_i, &m, &mm);
			_push_cigar1(&(map->cigar->cigar), &(map->cigar->cigar_n), &(map->cigar->cigar_m), (m+mm) << 4 | CMATCH);
			map->NM += mm;
		}
	}
    map->len_dif = bd - bi;
    map->bmax = bd > bi ? bd : bi;
	free(c); free(indel);
}

map_msg *map_init_msg(int n)
{
	int i, j, max_n = 1;
	map_msg *m_msg = (map_msg*)calloc(n, sizeof(map_msg));
    for (j = 0; j < n; ++j) {
        map_msg *p = m_msg+j;
        p->map_n = 0; p->map_m = max_n; p->map = (map_t*)malloc(max_n * sizeof(map_t));
        for (i = 0; i < max_n; ++i) {
            p->map[i].cigar = (cigar_t*)malloc(sizeof(cigar_t));
            p->map[i].cigar->cigar_n = 0; p->map[i].cigar->cigar_m = 10; p->map[i].cigar->cigar = (cigar32_t*)malloc(10 * sizeof(cigar32_t));
        }
    }
    return m_msg;
}

void map_free_msg(map_msg *m_msg, int n)
{
    int i, j;
    for (j = 0; j < n; ++j) {
        map_msg *p = m_msg + j;
        for (i = 0; i < p->map_m; ++i) {
            free(p->map[i].cigar->cigar);
            free(p->map[i].cigar);
        }
        free(p->map);
    }
    free(m_msg);
}

#ifdef _GEM_MAIN_
int main(int argc, char* argv[])
{
    char file[100]; FILE *mapf;
    strcpy(file, argv[1]);
    mapf = fopen(file, "r");
    map_msg *m_msg = map_init_msg(1); int max_n = 100;

	m_msg->map_n = 0;
	int line_len = LINE_SIZE, aln_len = LINE_SIZE;
	char *gem_line = (char*)malloc(LINE_SIZE * sizeof(char));
	char *aln_msg = (char*)malloc(LINE_SIZE * sizeof(char));
	char name[100]; int msg_len;
	char *t, *s; int t_i, _t_i;
	char chr[10], strand, md[200];
	long long offset; char os[100];
	int i;

	while (fgets(gem_line, line_len, mapf)) {
		sscanf(gem_line, "%[^\t]\t%*[^\t]\t%*[^\t]\t%[^\n]\n", name, aln_msg);
		if (strcmp(aln_msg, "-") == 0) {
			fprintf(stdout, "-\n");
			continue;
		}
		t = strtok(aln_msg, ",");
		msg_len = strlen(aln_msg);
		t_i = strlen(t)+1;
		while (t != NULL) {
			if (m_msg->map_n == m_msg->map_m) {
				m_msg->map_m <<= 1;
				m_msg->map = (map_t*)realloc(m_msg->map, m_msg->map_m * sizeof(map_t));
				if (m_msg->map == NULL) {
					fprintf(stderr, "[gem_map_read] Error: not enough memory.\n"); exit(0);
				}
				for (i = m_msg->map_m>>1; i < m_msg->map_m; ++i) {
					m_msg->map[i].cigar = (cigar_t*)malloc(sizeof(cigar_t));
					m_msg->map[i].cigar->cigar_n = 0; m_msg->map[i].cigar->cigar_m = 10;
					m_msg->map[i].cigar->cigar = (cigar32_t*)malloc(10 * sizeof(cigar32_t));
				}
			}
			_t_i = t_i;
			sscanf(t, "%[^:]:%c:%[^:]:%[^:]", chr, &strand, os, md);
            offset = atoll(os);
			t_i = _t_i;
			md2cigar(md, &(m_msg->map[m_msg->map_n]));
			m_msg->map[m_msg->map_n].offset = offset;
			m_msg->map[m_msg->map_n].strand = strand;
			strcpy(m_msg->map[m_msg->map_n].chr, chr);

			fprintf(stdout, "%s, %c, %lld, ", chr, strand, offset);
			printcigar(stdout, m_msg->map[m_msg->map_n].cigar->cigar, m_msg->map[m_msg->map_n].cigar->cigar_n);
			++(m_msg->map_n);
			if (t_i >= msg_len) break;
			t = strtok(aln_msg+t_i, ",");
			if (t == NULL) break;
			t_i += strlen(t)+1;
		}
	}

RET:
	if (gem_line) free(gem_line);
	if (aln_msg) free(aln_msg);

	map_free_msg(m_msg, 1);
	fclose(mapf);
	return 0;
}
#endif

int gem_map_read(FILE *mapf, map_msg *m_msg, char *gem_line, int line_size)
{
	if (fgets(gem_line, line_size, mapf) == NULL) return -1;
    
    int i, ct=0;
    gem_line[strlen(gem_line)-1] = 0;
    for (i = 0; i < strlen(gem_line); ++i) {
        if (gem_line[i] == '\t') {
            if (ct == 3) break;
            else ct++;
        }
    }
    m_msg->map_str = strdup(gem_line+i+1);
	return 0;
}

int gem_map_msg(map_msg *m_msg, int max_n)
{
	m_msg->map_n = 0;
    int i;
    char *aln_msg = m_msg->map_str; int msg_len;
	char *t; int t_i, _t_i;
	char chr[10], strand, md[200];
	long long offset; char os[100];

	//sscanf(gem_line, "%*[^\t]\t%*[^\t]\t%*[^\t]\t%[^\n]\n", aln_msg);
	if (strcmp(aln_msg, "-") == 0) {
#ifdef __DEBUG__
		fprintf(stderr, "-\n");
#endif
		goto RET;
	}

	msg_len = strlen(aln_msg);
	t = strtok(aln_msg, ",");
	t_i = strlen(t)+1;
	while (t != NULL) {
		if (m_msg->map_n >= max_n) {
			m_msg->map_n = 0;
			goto RET;
		}
		if (m_msg->map_n == m_msg->map_m) {
			m_msg->map_m <<= 1;
			m_msg->map = (map_t*)realloc(m_msg->map, m_msg->map_m * sizeof(map_t));
			if (m_msg->map == NULL) {
				fprintf(stderr, "[gem_map_read] Error: not enough memory.\n"); exit(0);
			}
			for (i = m_msg->map_m>>1; i < m_msg->map_m; ++i) {
				m_msg->map[i].cigar = (cigar_t*)malloc(sizeof(cigar_t));
				m_msg->map[i].cigar->cigar_n = 0; m_msg->map[i].cigar->cigar_m = 10;
				m_msg->map[i].cigar->cigar = (cigar32_t*)malloc(10 * sizeof(cigar32_t));
			}
		}
		_t_i = t_i;
		sscanf(t, "%[^:]:%c:%[^:]:%[^:]", chr, &strand, os, md);
        offset = atoll(os);
		t_i = _t_i;
		md2cigar(md, &(m_msg->map[m_msg->map_n]));
		if (strand == '-') _invert_cigar(&(m_msg->map[m_msg->map_n].cigar->cigar), m_msg->map[m_msg->map_n].cigar->cigar_n);
		m_msg->map[m_msg->map_n].offset = offset;
		m_msg->map[m_msg->map_n].strand = strand;
		strcpy(m_msg->map[m_msg->map_n].chr, chr);

#ifdef __DEBUG__
		fprintf(stderr, "%s, %c, %lld, ", chr, strand, offset);
		printcigar(stderr, m_msg->map[m_msg->map_n].cigar->cigar, m_msg->map[m_msg->map_n].cigar->cigar_n);
		fprintf(stderr, "\n");
#endif
		if (strcmp(chr, "chrM")!=0)
			++(m_msg->map_n);
		if (t_i >= msg_len) break;
		t = strtok(aln_msg+t_i, ",");
		if (t == NULL) break;
		t_i += strlen(t)+1;
	}
RET:
    free(m_msg->map_str);
	return m_msg->map_n;
}
