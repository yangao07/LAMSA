CC=			gcc
CFLAGS=		-g -Wall -O0
OBJS=		main.o build_ref.o bntseq.o lsat_aln.o frag_check.o split_mapping.o extend_ssw.o ssw.o ksw.o \
			./lsat_sam_parse/bam_aux.o ./lsat_sam_parse/bam.o ./lsat_sam_parse/bam_import.o ./lsat_sam_parse/kstring.o ./lsat_sam_parse/sam_header.o ./lsat_sam_parse/sam_view.o
PROG=		lsat
LIB=		-lm -lz
#MACRO=		-D __DEBUG__
MACRO=		-D SAM_IN
#MACRO=		-D PLAIN_IN
.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(MACRO) $< -o $@

all:$(PROG)
$(PROG):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIB) -o $@

clean:
	rm -f *.o lsat
