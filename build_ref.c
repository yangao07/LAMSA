/*
 * bulid_ref.c
 *
 * bulid index for reference sequence
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <getopt.h>
#include <zlib.h>
#include "bntseq.h"
//#include "bntseq_new.h"
#include "bwt.h"
#include "utils.h"
#include "build_ref.h"


int lsat_index_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   lsat index [option] <ref.fa>\n");
    fprintf(stderr, "                    bulid index for ref.fa\n\n");
    fprintf(stderr, "Option:  \n");
    fprintf(stderr, "         -i [INT]     Index program, <gem(0)>, <bwa(1)> or <soap2-dp(2)>. [Def=0]\n");
    fprintf(stderr, "\n");
    return 1;
}

void gem_build(char *prefix)
{
    char cmd[1024];
    sprintf(cmd, "./gem_index.sh %s", prefix);
    fprintf(stderr, "[lsat_index] Executing gem-indexer ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lsat_index] Indexing undone, gem-indexer exit abnormally.\n"); exit(1); }
    fprintf(stderr, "done.\n");
}

void bwa_build(char *prefix)
{
    char cmd[1024];
    sprintf(cmd, "./bwa_index.sh %s", prefix);
    fprintf(stderr, "[lsat_index] Executing bwa index ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lsat_index] Indexing undone, bwa index exit abnormally.\n"); exit(1); }
    fprintf(stderr, "done.\n");
}

void soap_bulid(char *prefix)
{
    char cmd[1024];
    sprintf(cmd, "./soap2dp_index.sh %s", prefix);
    fprintf(stderr, "[lsat_index] Executing soap2-dp-builder ... ");
    if (system(cmd) != 0 ) { fprintf(stderr, "\n[lsat_aln] Indexing undone, soap2-dp-builder exit abnormally.\n"); exit(1);}
    fprintf(stderr, " done.\n");
}

void bwt_index(char *prefix)
{
    int block_size = 10000000;
    int64_t l_pac;
    char *str, *str2, *str3;

	str  = (char*)calloc(strlen(prefix) + 10, 1);
	str2 = (char*)calloc(strlen(prefix) + 10, 1);
	str3 = (char*)calloc(strlen(prefix) + 10, 1);

    fprintf(stderr, "[bwt_index] Building bwt-index for genome...\n");
    { // for&rev.pac .ann .amb 
        gzFile fp = gzopen(prefix, "r");
        l_pac = bns_fasta2bntseq(fp, prefix, 0);
        gzclose(fp);
    }
    { // .bwt
        strcpy(str, prefix); strcat(str, ".pac");
        strcpy(str2, prefix); strcat(str2, ".bwt");
        bwt_bwtgen2(str, str2, block_size);

    }
    { // update .bwt
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".bwt");
        bwt = bwt_restore_bwt(str);
        bwt_bwtupdate_core(bwt);
        bwt_dump_bwt(str, bwt);
        bwt_destroy(bwt);
    }
    { // forward.pac
        gzFile fp = gzopen(prefix, "r");
        l_pac = bns_fasta2bntseq(fp, prefix, 1);
        gzclose(fp);
    }
    { // .sa
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".bwt");
        strcpy(str3, prefix); strcat(str3, ".sa");
        bwt = bwt_restore_bwt(str);
        bwt_cal_sa(bwt, 32);
        bwt_dump_sa(str3, bwt);
        bwt_destroy(bwt);
    }
    fprintf(stderr, "[bwt_index] Building done.\n");
    free(str); free(str2); free(str3);
}

int lsat_index(int argc, char *argv[])
{
	char *prefix=0;
	int type=0, c;
	//XXX clock_t t;  
	
	while ((c = getopt(argc, argv, "i:")) >= 0)
	{
		switch (c)
		{
            case 'i':
                type = atoi(optarg);
                if (type > 2 || type < 1) return lsat_index_usage();
                break;
            default:
                return lsat_index_usage();
		}
	}

	if (optind + 1 > argc) return lsat_index_usage();
	prefix = strdup(argv[optind]);

    bwt_index(prefix);
    if (type==0) gem_build(prefix);
    else if (type==1) bwa_build(prefix);
    else if(type==2) soap_bulid(prefix);
    else return lsat_index_usage();

	return 0;
}
