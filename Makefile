LIBS_PATH=libs

ifndef CC
	CC = gcc
endif

ifdef DEBUG
	OPT = -O0 -DDEBUG=1 --debug -g
else
	OPT = -O3
endif

CFLAGS = -Wall -Wextra $(OPT)

INCS=-I $(LIBS_PATH)/bit_array -I $(LIBS_PATH)/string_buffer \
     -I $(LIBS_PATH)/htslib/htslib -I $(LIBS_PATH)/seq_file -I src

LIBS=-L $(LIBS_PATH)/bit_array -L $(LIBS_PATH)/string_buffer \
     -L $(LIBS_PATH)/htslib/htslib -L src

LINK=-lalign -lstrbuf -lbitarr -lhts -lpthread -lz

# Compile and bundle all non-main files into library
CFILES=$(wildcard src/*.c)
ALIGN_FILES=$(filter-out src/*_cmdline.c,$(CFILES))
OBJ_FILES=$(ALIGN_FILES:.c=.o)

all: bin/needleman_wunsch bin/smith_waterman src/libalign.a examples

src/libalign.a: $(OBJ_FILES)
	ar -csru src/libalign.a $(OBJ_FILES)

%.o: %.c
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

bin/needleman_wunsch: bin src/nw_cmdline.c src/libalign.a
	$(CC) -o bin/needleman_wunsch $(CFLAGS) $(INCS) $(LIBS) src/nw_cmdline.c $(LINK)

bin/smith_waterman: bin src/sw_cmdline.c src/libalign.a
	$(CC) -o bin/smith_waterman $(CFLAGS) $(INCS) $(LIBS) src/sw_cmdline.c $(LINK)

bin:
	mkdir -p bin

examples: src/libalign.a
	cd examples; make

clean:
	rm -rf bin/* src/*.o src/libalign.a
	cd examples; make clean

.PHONY: all clean examples
