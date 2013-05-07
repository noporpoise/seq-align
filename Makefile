ifndef CC
	CC = gcc
endif

ifdef DEBUG
	OPT = -O0 -DDEBUG=1 --debug -g
else
	OPT = -O3
endif

CFLAGS = -Wall -Wextra $(OPT)

INCS = -I libs/bit_array -I libs/string_buffer \
       -I libs/htslib/htslib -I libs/seq_file/new_api -I src

LIBS = -lalign -lpthread -lz

LIB_OBJS=$(wildcard libs/bit_array/*.o) $(wildcard libs/string_buffer/*.o) \
         $(wildcard libs/htslib/htslib/*.o)

ALIGN_FILES=$(wildcard src/*.c)
OBJ_FILES=$(ALIGN_FILES:.c=.o)

all: bin/needleman_wunsch bin/smith_waterman src/libalign.a examples

src/libalign.a: $(OBJ_FILES)
	ar -csru src/libalign.a $(OBJ_FILES) $(LIB_OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

bin/needleman_wunsch: src/nw_cmdline.c src/libalign.a
	$(CC) -o bin/needleman_wunsch $(CFLAGS) $(INCS) -L src src/nw_cmdline.c $(LIBS)

bin/smith_waterman: src/sw_cmdline.c src/libalign.a
	$(CC) -o bin/smith_waterman $(CFLAGS) $(INCS) -L src src/sw_cmdline.c $(LIBS)

examples: src/libalign.a
	cd examples; make

clean:
	rm -rf bin/* src/*.o src/libalign.a
	cd examples; make clean

.PHONY: all clean examples
