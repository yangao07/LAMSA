# LAMSA
Long Approximated Matches-based Split Aligner

## Getting started
	git clone https://github.com/gaoyan07/LAMSA.git
	cd LAMSA; make
	./lamsa index ref.fa
	./lamsa aln ref.fa read.fq > aln.sam

## Introduction
LAMSA(Long Approximated Matches-based Split Aligner) is a  novel read alignment approach with faster speed and good ability of handling both co-linear and non-co-linear SV events.

LAMSA takes the advantage of the rareness of SVs to implement a specifically designed two-step split read alignment strategy, which efficiently solves the small events and mitigates the affection of repeats with co-linear alignment, and then well-handles the relatively large or non-co-linear events with a sparse dynamic programming (SDP)-based split alignment approach.

LAMSA has outstanding throughput on aligning both simulated and real datasets having various read length and sequencing error rates. It is severaly to over 100 folds faster than the state-of-art long read aligners. Morever, it also has good ability of handling various kinds of SV events within the read. 

LAMSA is open source and free for non-commercial use.

LAMSA is mainly designed by Bo Liu & Yan Gao and developed by Yan Gao in Center for Bioinformatics, Harbin Institute of Technology, China.

## Memory requirement
The memory usage of LAMSA can fit the configurations of most modern servers and workstations. Its peak memory footprint depends on the length of the read, i.e., 5.1 Gigabytes and 7.2 Gigabytes respectively for the 5000 bp and 100000 bp datasets, on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. These reads were aligned to GRCh37/hg19 reference genome.

## Installation
Current version of LAMSA needs to be run on Linux operating system.  
The source code is written in C, and can be directly download from: https://github.com/gaoyan07/LAMSA  
The makefile is attached. Use the make command for generating the executable file.  

## Synopsis

Index reference sequence and generate auxiliary files
```
lamsa index ref.fa
```
	
Align long read to reference
```
lamsa aln ref.fa read.fa/fq > aln.sam
```

## Commands and options
```
lamsa aln     [-t nThreads] [-l seedLen] [-i seedInv] [-p maxLoci] [-V maxSVLen]
              [-m matchScore] [-M mismatchScore] [-O gapOpenPen] [-E gapExtPen] 
              [-r maxOutputNum]  [-g minSplitLen] [-SC] [-o outSAM] 
               <ref.fa> <read.fa/fq>
              
Algorithm options:

    -t --thread    [INT]    Number of threads. [1]
    -l --seed-len  [INT]    Seed length. Moreover, LAMSA uses short sequence tool(e.g., GEM) to align
                            seeds and obtain their approximate matches. [50]
    -i --seed-inv  [INT]    Interval size of adjacent seeds. LAMSA extracts seeds on the starting
                            positons of every i bp. [100]
    -p --max-loci  [INT]    Maximum allowed number of a seed's locations. If a seed has more than -p
                            approximate matches, LAMSA would consider the seed is too repetitive, and
                            idiscard all the matches. [200]

    -V --SV-len    [INT]    Expected maximum length of SV. If the genomic distance of two seeds is
                            short than -V bp, they are avalibale to be connected to construct a
                            skeleton. [10000]

Scoring options:

    -m --match-sc  [INT]    Match score for SW-alignment. [1]
    -M --mis-pen   [INT]    Mismatch penalty for SW-alignment. [3]
    -O --open-pen  [INT]    Gap open penalty for SW-alignment. [5]
    -E --ext-pen   [INT]    Gap extension penalty for SW-alignment. A gap of length k costs O + k*E
                            (i.e. -O is for opening a zero-length gap). [2]

Output options:

    -r --max-out   [INT]    Maximum number of output records for a specific split read region. For a
                            specific region, LAMSA reserves the top -r alignment records. The record
                            with highest alignment score is considered as best alignment, others are
                            considered as alternative alignments. If the score of an alternative
                            alignment is less than half of the best alignment, it will not be output.
                            [10
    -g --gap-split [INT]    Minimum length of gap that causes a split-alignment. To avoid generating
                            insertion(I) or deletion(D) longer than -g bp in the SAM cigar. [100]

    -S --soft-clip          Use soft clipping for supplementary alignment. It is strongly recommended
                            to turn off this option to reduce the redundancy of output when mapping
                            relatively long reads. [false]
    -C --comment            Append FASTQ comment to SAM output. [false]
    -o --output    [STR]    Output file (SAM format). [stdout]

```

## Simulation benchmarking
We simulated a series of datasets from a variant human genomes. More precisely, we used RSVsim (version 1.10.0) to integrate 4002 SV events into human reference genome (GRCh37/hg19) to build the donor genome, including 532 duplications, 503 insertions, 2943 deletions and 24 inversions. The ratio of the four categories of SV events are configured by referring to the DGV database. Then, using the simulated donor genome as input, 15 datasets respectively with 5 kinds of read lengths (5000, 10000, 20000, 50000 and 100000 bp) and 3 kinds of sequencing error rates (1%, 2% and 4%) were simulated by wgsim(https://github.com/lh3/wgsim). Each of the datasets contains about 6 Giga bps, i.e., nearly 2X coverage of human genome. These datasets helped us to evaluate the performance of LAMSA. The donor genome and detailed list of simulated SV events have been uploaded to Google Drive, and can be downloaded through the following link: https://drive.google.com/folderview?id=0B24uQUND9m51UVlNWkJra19BMGs&usp=sharing


## Reference
LAMSA: fast split read alignment with long approximate matches. Manuscript in preparation.

## Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn or yangao07@hit.edu.cn


