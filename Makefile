CC      =	gcc
CFLAGS  =	-Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
DFLAGS  =	-g -Wall  
LIB     =	-lm -lz -lpthread

BIN_DIR =	.
SRC_DIR =   ./src

SOURCE  =	$(wildcard ${SRC_DIR}/*.c) 
#SOURCE  =	main.c build_ref.c bntseq.c lamsa_heap.c lamsa_aln.c lamsa_dp_con.c frag_check.c split_mapping.c ksw.c \
			gem_parse.c is.c bwtindex.c bwt_gen.c QSufSort.c kstring.c \
			bwt.c bwt_aln.c utils.c
OBJS    =	$(SOURCE:.c=.o)

BIN     =	$(BIN_DIR)/lamsa

GDB_DEBUG   =   $(BIN_DIR)/gdb_lamsa
NOR_DEBUG   =   $(BIN_DIR)/debug_lamsa
DMARCRO =	-D __DEBUG__

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

all:       $(SOURCE) $(BIN) 
#lamsa:     $(SOURCE) $(BIN) 
gdb_lamsa: $(SOURCE) $(GDB_DEBUG) 
debug_lamsa: $(SOURCE) $(NOR_DEBUG)


$(BIN): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LIB)

$(GDB_DEBUG):
	$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LIB)
$(NOR_DEBUG):
	$(CC) $(CFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LIB)

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)

clean_debug:
	rm -f $(SRC_DIR)/*.o $(GDB_DEBUG) $(NOR_DEBUG)
