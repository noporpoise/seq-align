LIBS_PATH=libs

ifdef DEBUG
	OPT = -O0 -DDEBUG=1 --debug -g -ggdb
else
	OPT = -O3
endif

CFLAGS = -Wall -Wextra -std=c99 $(OPT)
OBJFLAGS = -fPIC
LINKFLAGS = -lalign -lstrbuf -lbitarr -lpthread -lz

INCS=-I $(LIBS_PATH) -I src
LIBS=-L $(LIBS_PATH)/bit_array -L $(LIBS_PATH)/string_buffer -L src
LINK=-lalign -lstrbuf -lbitarr -lpthread -lz

# If we are missing libraries, force (cd libs && make)
REQ=$(LIBS_PATH)/bit_array/Makefile $(LIBS_PATH)/string_buffer/Makefile \
    $(LIBS_PATH)/seq_file/Makefile $(LIBS_PATH)/sort_r/Makefile

# Compile and bundle all non-main files into library
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)

all: bin/needleman_wunsch bin/smith_waterman bin/lcs src/libalign.a examples

src/libalign.a: $(OBJS)
	ar -csru src/libalign.a $(OBJS)

$(LIBS_PATH)/%/Makefile:
	cd libs; make;

$(OBJS): $(REQ)

%.o: %.c
	$(CC) $(CFLAGS) $(OBJFLAGS) $(INCS) -c $< -o $@

bin/needleman_wunsch: src/tools/nw_cmdline.c src/libalign.a | bin
	$(CC) -o bin/needleman_wunsch $(SRCS) $(TGTFLAGS) $(INCS) $(LIBS) src/tools/nw_cmdline.c $(LINKFLAGS)

bin/smith_waterman: src/tools/sw_cmdline.c src/libalign.a | bin
	$(CC) -o bin/smith_waterman $(SRCS) $(TGTFLAGS) $(INCS) $(LIBS) src/tools/sw_cmdline.c $(LINKFLAGS)

bin/lcs: src/tools/lcs_cmdline.c src/libalign.a | bin
	$(CC) -o bin/lcs $(SRCS) $(TGTFLAGS) $(INCS) $(LIBS) src/tools/lcs_cmdline.c $(LINKFLAGS)

bin/seq_align_tests: src/tools/tests.c src/libalign.a
	mkdir -p bin
	$(CC) -o $@ $< $(CFLAGS) $(INCS) $(LIBS) $(LINK)

examples: src/libalign.a
	cd examples; make LIBS_PATH=$(abspath $(LIBS_PATH))

bin:
	mkdir -p bin

clean:
	rm -rf bin src/*.o src/libalign.a
	cd examples && make clean

test: bin/seq_align_tests
	./bin/seq_align_tests

.PHONY: all clean examples test
