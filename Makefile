CC=			gcc
CFLAGS=		-g -Wall -O0
OBJS=		main.o build_ref.o	bntseq.o lsat_aln.o
PROG=		lsat
LIB=		-lm -lz
.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

all:$(PROG)
lsat:$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIB) -o $@

clean:
	rm -f *.o lsat
