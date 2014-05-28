CC ?= gcc
LIBS_PATH=libs

PLATFORM := $(shell uname)
COMPILER := $(shell ($(CC) -v 2>&1) | tr A-Z a-z )

ifdef DEBUG
	OPT = -O0 -DDEBUG=1 --debug -g -ggdb
else
	ifneq (,$(findstring clang,$(COMPILER)))
		# clang Link Time Optimisation (lto) seems to have issues atm
		OPT = -O3
	else
		OPT = -O4
		#OPT = -O4 -flto
		#TGTFLAGS = -fwhole-program
	endif
endif

CFLAGS = -Wall -Wextra $(OPT)
OBJFLAGS = -fPIC
LINKFLAGS = -lalign -lstrbuf -lbitarr -lpthread -lz

INCS=-I $(LIBS_PATH)/bit_array -I $(LIBS_PATH)/string_buffer \
     -I $(LIBS_PATH)/seq_file -I src

LIBS=-L $(LIBS_PATH)/bit_array -L $(LIBS_PATH)/string_buffer -L src

REQ=$(LIBS_PATH)/bit_array/Makefile $(LIBS_PATH)/string_buffer/Makefile $(LIBS_PATH)/seq_file/Makefile

# Compile and bundle all non-main files into library
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)

all: bin/needleman_wunsch bin/smith_waterman src/libalign.a examples

src/libalign.a: $(OBJS)
	ar -csru src/libalign.a $(OBJS)

$(LIBS_PATH)/%/Makefile:
	cd libs; make;

$(OBJS): $(REQ)

%.o: %.c
	$(CC) $(CFLAGS) $(OBJFLAGS) $(INCS) -c $< -o $@

bin/needleman_wunsch: bin tools/nw_cmdline.c src/libalign.a
	$(CC) -o bin/needleman_wunsch $(SRCS) $(TGTFLAGS) $(INCS) $(LIBS) tools/nw_cmdline.c $(LINKFLAGS)

bin/smith_waterman: bin tools/sw_cmdline.c src/libalign.a
	$(CC) -o bin/smith_waterman $(SRCS) $(TGTFLAGS) $(INCS) $(LIBS) tools/sw_cmdline.c $(LINKFLAGS)

bin:
	mkdir -p bin

examples: src/libalign.a
	cd examples; make LIBS_PATH=$(abspath $(LIBS_PATH))

clean:
	rm -rf bin src/*.o src/libalign.a
	cd examples; make clean

.PHONY: all clean examples
