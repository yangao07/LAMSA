CC=			gcc
CFLAGS=		-g -Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
OBJS=		main.o build_ref.o bntseq.o lsat_heap.o lsat_aln.o lsat_dp_con.o frag_check.o split_mapping.o ksw.o \
			gem_parse.o is.o bwtindex.o bwt_gen.o QSufSort.o\
			./lsat_sam_parse/bam_aux.o ./lsat_sam_parse/bam.o ./lsat_sam_parse/bam_import.o \
			./lsat_sam_parse/kstring.o ./lsat_sam_parse/sam_header.o ./lsat_sam_parse/sam_view.o \
			bwt.o bwt_aln.o utils.o
#PROG=		lsat_gdb
PROG=		lsat
PROG1=      ~/bin/lsat
LIB=		-lm -lz -lpthread
#MACRO=     -D __NEW__
#MACRO=		-D __DEBUG__
.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(MACRO) $< -o $@

all:$(PROG) $(PROG1)
$(PROG):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIB) -o $@
$(PROG1):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIB) -o $@

clean:
	rm -f *.o lsat lsat_gdb ~/bin/lsat
