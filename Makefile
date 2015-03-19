CC=			gcc
CFLAGS=		-g -Wall -O0 -Wno-unused-variable -Wno-unused-but-set-variable
OBJS=		main.o build_ref.o bntseq.o lsat_heap.o lsat_aln.o frag_check.o split_mapping.o extend_ssw.o ssw.o ksw.o \
			gem_parse.o \
			./lsat_sam_parse/bam_aux.o ./lsat_sam_parse/bam.o ./lsat_sam_parse/bam_import.o \
			./lsat_sam_parse/kstring.o ./lsat_sam_parse/sam_header.o ./lsat_sam_parse/sam_view.o \
			bwt.o bwt_aln.o utils.o
PROG=		lsat
LIB=		-lm -lz
#MACRO=		-D __DEBUG__
.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(MACRO) $< -o $@

all:$(PROG)
$(PROG):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIB) -o $@

clean:
	rm -f *.o lsat
