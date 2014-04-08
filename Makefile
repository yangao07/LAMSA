CC=			gcc
CFLAGS=		-g -Wall -O0
OBJS=		main.o build_ref.o bntseq.o lsat_aln.o frag_check.o split_mapping.o extend_ssw.o ssw.o ksw.o
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
