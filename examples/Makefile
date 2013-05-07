CFLAGS=-Wall -Wextra

all: nw_example sw_example

sw_example: sw_example.c ../libalign.a
	$(CC) $(CFLAGS) -I ../src -I ../libs/bit_array/ -L .. -L ../libs/bit_array/ \
	sw_example.c -o sw_example -lalign -lbitarr

nw_example: nw_example.c ../libalign.a
	$(CC) $(CFLAGS) -I ../src -L .. nw_example.c -o nw_example -lalign

clean:
	rm -rf nw_example sw_example *.greg *.dSYM

.PHONY: all clean
