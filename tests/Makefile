CC ?= gcc
LIBS_PATH=../libs
CFLAGS=-Wall -Wextra -std=c99

INCS=-I $(LIBS_PATH)/bit_array -I $(LIBS_PATH)/string_buffer \
     -I $(LIBS_PATH)/seq_file -I ../src

LIBS=-L $(LIBS_PATH)/bit_array -L $(LIBS_PATH)/string_buffer -L ../src

LINK=-lalign -lstrbuf -lbitarr -lpthread -lz

all: tests
	
tests: tests.c ../src/libalign.a
	$(CC) tests.c -o tests $(CFLAGS) $(INCS) $(LIBS) $(LINK)

clean:
	rm -rf tests *.greg *.dSYM

.PHONY: all clean