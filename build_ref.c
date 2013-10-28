/*
 * bulid_ref.c
 *
 * bulid index for reference sequence(soap2-dp)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <zlib.h>
#include "bntseq.h"
#include "build_ref.h"

int soap_bulid(char *prefix, char *path)
{
	if (prefix != 0 && path != NULL)
	{
		//soap2-dp-bulid prefix
		char cmd[1024];
		
		strcpy(cmd, path); strcat(cmd, "/soap2-dp-builder "); strcat(cmd, prefix);
        strcat(cmd, " > "); strcat(cmd, prefix), strcat(cmd, ".build");
		fprintf(stderr, "Executing soap2-dp-builder ... \n");
		fprintf(stderr, "%s", cmd);
		if (system(cmd) != 0 )
			exit(-1);
		fprintf(stderr, " done.\n");
		return 1;
	}
	return 0;
}
int lsat_index(int argc, char *argv[])
{
	char *prefix=0;
	char path[1024]="\0";
	int c;
	//XXX clock_t t;  
	
	while ((c = getopt(argc, argv, "s:")) >= 0)
	{
		switch (c)
		{
			case 's':
				strcpy(path, optarg);
				break;
			default:
				return 1;
		}
	}

	if (optind + 1 > argc)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   lsat index [-s] <ref.fa>\n");
		fprintf(stderr, "                    bulid soap2-dp index for ref.fa\n\n");
		fprintf(stderr, "Option:  -s STR     soap2-dp's path. [Def=\"./soap2-dp/\"]\n\n");
		return 1;
	}
	prefix = strdup (argv[optind]);
	if (strlen(path) == 0)
	{
		if (getcwd(path, sizeof(path)) == NULL)
		{
			perror("getcwd error");
			exit(-1);
		}
		strcat(path, "/soap2-dp");
	}

	//soap2-dp-builder
	soap_bulid(prefix, path);

	//lsat_fa2pac
	gzFile fp = gzopen(prefix, "r");
	bns_fa2bnt(fp, prefix);

	free(prefix);
	return 0;
}
