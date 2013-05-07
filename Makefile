LIBS_PATH=./libs

STRING_BUF_PATH = $(LIBS_PATH)/string_buffer/
BIT_ARRAY_PATH = $(LIBS_PATH)/bit_array/
SEQ_FILE_PATH = $(LIBS_PATH)/seq_file/new_api/
HTS_PATH = $(LIBS_PATH)/htslib/

ifndef CC
	CC = gcc
endif

ifdef DEBUG
	CFLAGS := -O0 -DDEBUG=1 --debug -g
else
	CFLAGS := -O3
endif

override HTS_PATH:=$(HTS_PATH)/htslib/

# Add data type for alignment scoring
CFLAGS := $(CFLAGS) -Wall -Wextra

LIBPATHS = -L $(BIT_ARRAY_PATH) -L $(STRING_BUF_PATH) -L $(HTS_PATH) -L .

INCS = -I $(BIT_ARRAY_PATH) -I $(STRING_BUF_PATH) \
       -I $(HTS_PATH) -I $(SEQ_FILE_PATH) -I src

LIBS = -lstrbuf -lbitarr -lhts -lpthread -lz

ALIGN_FILES=$(wildcard src/*.c) src/sort_r.c src/needleman_wunsch.c src/smith_waterman.c
OBJ_FILES=$(ALIGN_FILES:.c=.o)

all: bin/needleman_wunsch bin/smith_waterman libalign.a examples

libalign.a: $(OBJ_FILES)
	ar -csru libalign.a $(OBJ_FILES)

%.o: %.c
	$(CC) $(CFLAGS) $(OPT) $(INCS) -c $< -o $@

bin/needleman_wunsch: src/nw_cmdline.c libalign.a
	$(CC) -o bin/needleman_wunsch $(CFLAGS) $(INCS) $(LIBPATHS) src/nw_cmdline.c -lalign $(LIBS)

bin/smith_waterman: src/sw_cmdline.c libalign.a
	$(CC) -o bin/smith_waterman $(CFLAGS) $(INCS) $(LIBPATHS) src/sw_cmdline.c -lalign $(LIBS)

examples: libalign.a
	cd examples; make

clean:
	rm -rf bin/* src/*.o libalign.a
	cd examples; make clean

.PHONY: all clean examples
