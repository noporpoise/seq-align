CC ?= gcc
LIBS_PATH=../libs
CFLAGS=-Wall -Wextra -std=c99

INCS=-I $(LIBS_PATH)/bit_array -I $(LIBS_PATH)/string_buffer \
     -I $(LIBS_PATH)/seq_file -I ../src

LIBS=-L $(LIBS_PATH)/bit_array -L $(LIBS_PATH)/string_buffer -L ../src

LINK=-lalign -lstrbuf -lbitarr -lpthread -lz

all: nw_example sw_example lcs_example

sw_example: sw_example.c ../src/libalign.a
	$(CC) sw_example.c -o sw_example $(CFLAGS) $(INCS) $(LIBS) $(LINK)

nw_example: nw_example.c ../src/libalign.a
	$(CC) nw_example.c -o nw_example $(CFLAGS) $(INCS) $(LIBS) $(LINK)
	
lcs_example: lcs_example.c ../src/libalign.a
	$(CC) lcs_example.c -o lcs_example $(CFLAGS) $(INCS) $(LIBS) $(LINK)

clean:
	rm -rf nw_example sw_example lcs_example *.greg *.dSYM

.PHONY: all clean
