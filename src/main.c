#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "lamsa_aln.h"
#include "lamsa_index.h"

#define VERSION "1.0.0"
#define CONTACT "Yan Gao <yangao07@hit.edu.cn>"

static int usage(void)	//main usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: LAMSA (Long Approximated Matches-based Split Aligner)\n");
    fprintf(stderr, "Version: %s\n", VERSION);
    fprintf(stderr, "Contact: %s\n\n", CONTACT); 

	fprintf(stderr, "Usage:   lamsa <command> [options]\n\n");

	fprintf(stderr, "Command: \n");
	fprintf(stderr, "         index       index reference sequence\n");
	fprintf(stderr, "         aln         align long sequence to reference\n");

	fprintf(stderr, "\n");
	return 1;
}

char lamsa_pg[1024];
int main(int argc, char *argv[])
{
	int i;
    sprintf(lamsa_pg, "@PG\tID:lamsa\tPN:lamsa\tVN:%s\tCL:%s", VERSION, argv[0]);
	for (i = 1; i < argc; ++i) sprintf(lamsa_pg+strlen(lamsa_pg), " %s", argv[i]);
    
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0)      return lamsa_index(argc, argv);
	else if (strcmp(argv[1], "aln") == 0)   return lamsa_aln(argc, argv);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}

    free(lamsa_pg);
	return 0;
}
