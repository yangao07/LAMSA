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

static int usage(void)	//main usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: lsat (A Long Sequence Alignment Tool.)\n");
	fprintf(stderr, "Usage:   lsat <cmd> [opts]\n\n");
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "         index      index reference sequence(soap2-dp)\n");
	fprintf(stderr, "         aln        align long sequence to reference\n");

	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0)      return lsat_index(argc-1, argv+1);
	else if (strcmp(argv[1], "aln") == 0)   return lsat_aln(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
