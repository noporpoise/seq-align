LIBS_PATH=../libs

ifndef CC
	CC = gcc
endif

ifdef DEBUG
	CFLAGS := -DDEBUG=1 --debug -g
else
	CFLAGS := -O3
endif

UTILITY_LIB_PATH = $(LIBS_PATH)/utility_lib
STRING_BUF_PATH = $(LIBS_PATH)/string_buffer
BIT_ARRAY_PATH = $(LIBS_PATH)/bit_array
SEQ_FILE_PATH = $(LIBS_PATH)/seq_file
SAMTOOLS_PATH = $(HOME)/bioinf/samtools-0.1.18

# Add data type for alignment scoring
CFLAGS := $(CFLAGS) -Wall -Wextra \
          -I. -Icommon/ -I$(SEQ_FILE_PATH) -I$(UTILITY_LIB_PATH) \
          -I$(BIT_ARRAY_PATH) -I$(STRING_BUF_PATH) -I$(SAMTOOLS_PATH)

LIB_INCS = -L$(SEQ_FILE_PATH) -L$(UTILITY_LIB_PATH) \
           -L$(BIT_ARRAY_PATH) -L$(STRING_BUF_PATH) -L$(SAMTOOLS_PATH) -L.

LIB_LIST = -lseqfile -lstrbuf -lbitarr -lbam -lutil -lz

NW_ARGS = -DSCORE_TYPE='int'
SW_ARGS = -DSCORE_TYPE='unsigned int'

NW_FILES = nw_cmdline.c needleman_wunsch.c mem_size.c common/*.c
SW_FILES = sw_cmdline.c smith_waterman.c mem_size.c common/*.c

all: clean
	$(CC) -o needleman_wunsch $(CFLAGS) $(LIB_INCS) $(NW_ARGS) $(NW_FILES) $(LIB_LIST)
	$(CC) -o smith_waterman $(CFLAGS) $(LIB_INCS) $(SW_ARGS) $(SW_FILES) $(LIB_LIST)

clean:
	rm -rf needleman_wunsch needleman_wunsch.dSYM needleman_wunsch.greg
	rm -rf smith_waterman smith_waterman.dSYM smith_waterman.greg
