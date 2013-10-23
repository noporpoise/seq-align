LIBS_PATH=../libs
CFLAGS=-Wall -Wextra

INCS=-I $(LIBS_PATH)/bit_array -I $(LIBS_PATH)/string_buffer \
     -I $(LIBS_PATH)/seq_file -I ../src

LIBS=-L $(LIBS_PATH)/bit_array -L $(LIBS_PATH)/string_buffer -L ../src

LINK=-lalign -lstrbuf -lbitarr -lpthread -lz

all: nw_example sw_example

sw_example: sw_example.c ../src/libalign.a
	$(CC) sw_example.c -o sw_example $(CFLAGS) $(INCS) $(LIBS) $(LINK)

nw_example: nw_example.c ../src/libalign.a
	$(CC) nw_example.c -o nw_example $(CFLAGS) $(INCS) $(LIBS) $(LINK)

../src/libalign.a:
	cd ..; make src/libalign.a

clean:
	rm -rf nw_example sw_example *.greg *.dSYM

.PHONY: all clean
