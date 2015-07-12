/*
 * main.c
 *
 * LSAT: A Long Sequence Alignment Tool.
 *
 * Created by Yan Gao on 2013/9/20.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "lsat_aln.h"
#include "build_ref.h"

#define VERSION "1.0.0"

static int usage(void)	//main usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: lsat (A Long Sequence Alignment Tool.)\n");
	fprintf(stderr, "Usage:   lsat <cmd> [opts]\n\n");
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "         index      index reference sequence\n");
	fprintf(stderr, "         aln        align long sequence to reference\n");

	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
    extern char *lsat_pg;
    char *pg = (char*)malloc(1024*sizeof(char));
	int i;
    sprintf(pg, "@PG\tID:lsat\tPN:lsat\tVN:%s\tCL:%s", VERSION, argv[0]);
	for (i = 1; i < argc; ++i) sprintf(pg+strlen(pg), " %s", argv[i]);
    lsat_pg = pg;
    
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0)      return lsat_index(argc-1, argv+1);
	else if (strcmp(argv[1], "aln") == 0)   return lsat_aln(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}

    free(lsat_pg);
	return 0;
}
