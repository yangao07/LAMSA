/*
 * extend_main.c
 *
 * Construct ref and read seq of fragment.
 * Then call extend_ssw().
 *
 * Created by Yan Gao on 08/26/13.
 */
//contig 1620 read_n 494
//1       + 247142975 - 810000 500
//		    247115374 - 809862 362
//1       + 247115200 - 809861 361
//		    247064398 - 809607 107
//1       + 247064444 - 809606 106
//		    247064444 - 809606 106
//1       + 247064203 - 809605 105
//	        247064203 - 809605 105
//1       + 247064085 - 809604 104
//		    247043485 - 809501 1

/*
 * ref:  ./data/chr1_hg18.fa
 * read: ./data/venter_chr1_100000_0.01_het.fq 
 */
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
//#include "ssw.h"
//#include "kseq.h"

//KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{

	frag_check("../data/frag.msg");
	return 0;
}

