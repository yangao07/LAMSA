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
#include "lamsa_index.h"


int lamsa_index_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   lamsa index <ref.fa>\n");
    fprintf(stderr, "                     bulid BWT and GEM index for ref.fa\n\n");
    //fprintf(stderr, "Option:  \n");
    //fprintf(stderr, "         -i [INT]     Index program, <gem(0)>, <bwa(1)> or <soap2-dp(2)>. [Def=0]\n");
    fprintf(stderr, "\n");
    return 1;
}

char *get_bin_dir(char *bin)
{
    char *end = strrchr(bin, '/');
    bin[end-bin] = '\0';
    return bin;
}

void gem_build(char *bin_dir, char *prefix)
{
    char cmd[1024];

    sprintf(cmd, "bash %s/gem/gem_index.sh %s", bin_dir, prefix);
    //sprintf(cmd, "./gem/gem-indexer -i %s -o %s 2> /dev/null", prefix, prefix);
    fprintf(stderr, "[lamsa_index] Executing gem-indexer ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lamsa_index] Indexing undone, gem-indexer exit abnormally.\n"); exit(1); }
    fprintf(stderr, "done!\n");
}

void bwa_build(char *bin_dir, char *prefix)
{
    char cmd[1024];
    sprintf(cmd, "bash %s/bwa/bwa_index.sh %s", bin_dir, prefix);
    fprintf(stderr, "[lamsa_index] Executing bwa index ... ");
    if (system(cmd) != 0) { fprintf(stderr, "\n[lamsa_index] Indexing undone, bwa index exit abnormally.\n"); exit(1); }
    fprintf(stderr, "done!\n");
}

void soap_bulid(char *bin_dir, char *prefix)
{
    char cmd[1024];
    sprintf(cmd, "bash %s/soap2-dp/soap2dp_index.sh %s", bin_dir, prefix);
    fprintf(stderr, "[lamsa_index] Executing soap2-dp-builder ... ");
    if (system(cmd) != 0 ) { fprintf(stderr, "\n[lamsa_aln] Indexing undone, soap2-dp-builder exit abnormally.\n"); exit(1);}
    fprintf(stderr, " done!\n");
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
    fprintf(stderr, "[bwt_index] Building done!\n");
    free(str); free(str2); free(str3);
}

int lamsa_index(int argc, char *argv[])
{
    clock_t t = clock();
	char *prefix=0;
	int type=0, c;
	
	/*while ((c = getopt(argc, argv, "i:")) >= 0)
	{
		switch (c)
		{
            case 'i':
                type = atoi(optarg);
                if (type > 2 || type < 1) return lamsa_index_usage();
                break;
            default:
                return lamsa_index_usage();
		}
	}*/

	if (optind + 2 > argc) return lamsa_index_usage();
	prefix = strdup(argv[optind+1]);

    char *bin_dir;
    bin_dir = get_bin_dir(argv[0]);

    bwt_index(prefix);
    if (type==0) gem_build(bin_dir, prefix);
    else if (type==1) bwa_build(bin_dir, prefix);
    else if(type==2) soap_bulid(bin_dir, prefix);
    else return lamsa_index_usage();
    
    //fprintf(stderr, "[lamsa_index] Time Consupmtion: %.3f sec.\n", (float)(clock()-t)/CLOCKS_PER_SEC);
	return 0;
}
