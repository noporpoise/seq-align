CFLAGS=-Wall -Wextra

all: nw_example sw_example

sw_example: sw_example.c ../src/libalign.a
	$(CC) $(CFLAGS) -I ../src -L ../src sw_example.c -o sw_example -lalign

nw_example: nw_example.c ../src/libalign.a
	$(CC) $(CFLAGS) -I ../src -L ../src nw_example.c -o nw_example -lalign

../src/libalign.a:
	cd ..; make src/libalign.a

clean:
	rm -rf nw_example sw_example *.greg *.dSYM

.PHONY: all clean
