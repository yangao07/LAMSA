/*
 * main.c
 *
 * LSAT: A Long Sequence Alignment Tool.
 *
 * Created by Yan Gao on 2013/9/20
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

static struct option long_options[] = {
    { "help", 0, NULL, 0 },
    { "ref-prefix",  1, NULL, 'r' },
    { "no-soap2-dp", 0, NULL, 'n' },
    { "aln-result",  1, NULL, 'a' },
    { "max-error",   1, NULL, 'm' },
    { "per-len",     1, NULL, 'l' },
    { "frag-out",    1, NULL, 'f' },
    { "ctg-msg",     1, NULL, 'c' },
	{ "output",      1, NULL, 'o' },
    { 0,             0,    0,   0 },
};
/*
static int usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: %s (A Long Sequence Alignment Tool.)\n", argv[0]);
    fprintf(stderr, "Usage:   %s <soap2-dp index> <read file> [options]\n\n", argv[0]);
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "         --help    Print this usage.\n");
    //fprintf(stderr, "         -r --ref-prefix [STR]  The prefix of soap2-dp reference index.\n");
    fprintf(stderr, "         -n --no-soap2-dp\n");
    fprintf(stderr, "                   With soap2-dp alignment result existed already.\n");
    fprintf(stderr, "         -a --aln-reslut [STR]\n");
    fprintf(stderr, "                   The soap2-dp alignment result. Only if '-n' is used, [Def=\"prefix.out.0\"]\n");
    fprintf(stderr, "         -m --max-error[INT][STR] (soap2-dp option)\n");
    fprintf(stderr, "                   Soap2-dp options. Maximum #errors allowed. [Def=3e]\n");
    //fprintf(stderr, "         -h [INT] Alignment type. [Def=4 Random Best Alignment]\n");
    fprintf(stderr, "         -l --per-len [INT]\n");
    fprintf(stderr, "                   The length of per-read. [Def=100]\n");
    //fprintf(stderr, "         -f --frag-out [STR]\n");
    //fprintf(stderr, "                   Output the fragment message to file. [Def=\"prefix_frag.msg\"]\n");
    //fprintf(stderr, "         -c --ctg-msg [STR]\n");
    //fprintf(stderr, "                   Output the congit message to file. [Def=\"prefix_ctg.msg\"]\n");
	fprintf(stderr, "         -o --output [STR]\n");
	fprintf(stderr, "                     Output the alignment result to file(SAM format). [Def=\"prefix_out.sam\"]\n");

	return 0;
}
*/
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
