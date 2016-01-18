# LAMSA
Long Approximate Matches-based Split Aligner

## Getting started
	git clone https://github.com/hitbc/LAMSA.git
	cd LAMSA; make
	./lamsa index ref.fa
	./lamsa aln ref.fa read.fq > aln.sam

## Introduction
LAMSA (Long Approximate Matches-based Split Aligner) is a  novel split alignment approach with faster speed and good ability of handling SV events. It is well-suited to align long reads (over thousands of base-pairs). 

LAMSA takes the advantage of the rareness of SVs to implement a specifically designed two-step split read alignment strategy. That is, LAMSA initially splits the read into fragments and co-linearly align the fragments to solve the small events and mitigate the affection of repeats. The alignments of the fragments are then used for implementing a sparse dynamic programming (SDP)-based split alignment approach to further handle the large or non-co-linear events.

LAMSA has outstanding throughput on aligning both simulated and real datasets having various read length and sequencing error rates. It is several to over 100 folds faster than the state-of-art long read aligners. Moreover, it also has good ability of handling various kinds of SV events within the read. 

LAMSA is open source and free for non-commercial use.

LAMSA is mainly designed by Bo Liu & Yan Gao and developed by Yan Gao in Center for Bioinformatics, Harbin Institute of Technology, China.

## Memory usage
The memory usage of LAMSA can fit the configurations of most modern PCs. Its peak memory footprint depends on the length of the read, i.e., 5.1 Gigabytes and 7.2 Gigabytes respectively for the 5000 bp and 100000 bp datasets, on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. These reads were aligned to human reference genome GRCh37/hg19.

## Installation
Current version of LAMSA needs to be run on Linux operating system.  
The source code is written in C, and can be directly download from: https://github.com/hitbc/LAMSA
A mirror is also in: https://github.com/gaoyan07/LAMSA  
Moreover, in current version of LAMSA, we employed the GEM mapper (http://gemlibrary.sourceforge.net/) for generating the approximate matches of the fragments of reads. To be more user-friendly, we have built the source code of GEM mapper (version core_i3-20130406-045632) into that of LAMSA. In both of the genome indexing and read alignment, LAMSA will call the corresponding functions of GEM mapper automatically.
The makefile of LAMSA is attached. Use the make command for generating the executable file.  

## Synopsis

Reference genome indexing
```
lamsa index ref.fa
```
	
Read alignment
```
lamsa aln ref.fa read.fa/fq > aln.sam
```

## Commands and options
```
lamsa aln     [-t nThreads] [-l seedLen] [-i seedInv] [-p maxLoci] [-V maxSVLen] 
              [-v overlapRatio] [-s maxSkeletonNum] [-R bwtMaxReg] [-k bwtKmerLen]
              [-m matchScore] [-M mismatchScore] [-O gapOpenPen] [-E gapExtPen] 
              [-r maxOutputNum]  [-g minSplitLen] [-SC] [-o outSAM] 
              <ref.fa> <read.fa/fq>
              
Algorithm options:

    -t --thread    [INT]    Number of threads. [1]
    -l --seed-len  [INT]    Length of seeding fragments. Moreover, LAMSA splits the read into a 
                            series of -l bp long fragments, and employs NGS aligner to generate the 
                            approximate matches of the fragments. [50]
    -i --seed-inv  [INT]    Distance between neighboring seeding fragments. LAMSA extracts seeding 
                            fragments starting at every -i bp of the read. [100]
    -p --max-loci  [INT]    Maximum allowed number of hits. If a seeding fragment has more than -p 
                            approximate matches, LAMSA would consider the seed is too repetitive, and
                            discard all the matches. [200]
    -V --SV-len    [INT]    Expected maximum length of SV. If the genomic distance of two seeding 
                            fragments is longer than -V bp, they cannot be connected to build a legal 
                            edge in the sparse dynamic programming process. [10000]
    -v --ovlp-rat  [FLOAT]  Minimum overlapping ratio to cluster two skeletons or alignment records. 
                            [0.7]
    -s --max-skel  [INT]    Maximum number of skeletons that are reserved in a cluster for a specific 
                            read region. For a specific region of read, LAMSA reserves the top -s 
                            skeletons. These skeletons are used to generate best and alternative 
                            alignment records. [10]
    -R --max-reg   [INT]    Maximum allowed length of unaligned read part to trigger a bwt-based query. 
                            For a read part being unaligned after all the sparse dynamic programming 
                            process-based split alignment, if it is longer than -R bp, LAMSA will not 
                            further process it; otherwise, LAMSA will query the exact matches of the 
                            k-mers of the unaligned part as hits to further align the read part. [300]
    -k --bwt-kmer  [INT]    Length of BWT-seed. For the unaligned read part shorter than -R bp, LAMSA 
                            will extracts all its -k bp tokens and query their exact matches as hits. 
                            [19]

Scoring options:

    -m --match-sc  [INT]    Match score for SW-alignment. [1]
    -M --mis-pen   [INT]    Mismatch penalty for SW-alignment. [3]
    -O --open-pen  [INT]    Gap open penalty for SW-alignment. [5]
    -E --ext-pen   [INT]    Gap extension penalty for SW-alignment. A gap of length k costs O + k*E.
                            [2]

Output options:
    -r --max-out   [INT]    Maximum number of output records for a specific split read region. For a 
                            specific region, LAMSA reserves the top -r alignment records. The record 
                            with highest alignment score is considered as best alignment, others are 
                            considered as alternative alignments. Moreover, if the score of an 
                            alternative alignment is less than half of the best alignment, it will not 
                            be output. [10]
    -g --gap-split [INT]    Minimum length of gap that causes a split-alignment. To avoid generating
                            insertion(I) or deletion(D) longer than -g bp in the SAM cigar. [100]
    -S --soft-clip          Use soft clipping for supplementary alignment. It is strongly recommended
                            to turn off this option to reduce the redundancy of output when mapping
                            relatively long reads. [false]
    -C --comment            Append FASTQ comment to SAM output. [false]
    -o --output    [STR]    Output file (SAM format). [stdout]

```

## Simulation benchmarking
We simulated a series of datasets from an in silico donor human genome. More precisely, we used RSVsim (version 1.10.0) to integrate 4002 SV events into human reference genome (GRCh37/hg19) to build the donor genome, including 532 duplications, 503 insertions, 2943 deletions and 24 inversions. The maximum size of the SV events are 10000 bp, and minimum size of the duplications, deletions and insertions is 50 bp, while the minimum size of the inversions is 500 bp. The ratio and the size of the SV events are configured by referring to the DGV database (http://dgv.tcag.ca/dgv/app/home). Then, using the simulated donor genome as input, 15 datasets respectively with 5 kinds of read lengths (5000, 10000, 20000, 50000 and 100000 bp) and 3 kinds of sequencing error rates (1%, 2% and 4%) were simulated by Wgsim(https://github.com/lh3/wgsim). Each of the datasets contains about 6 Giga bps. These datasets helped us to evaluate the performance of LAMSA. The donor genome and detailed list of simulated SV events are available at: https://sourceforge.net/projects/lamsa/files/data/


## Reference
LAMSA: fast split read alignment with long approximate matches. Manuscript submitted.

## Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn or yangao07@hit.edu.cn


