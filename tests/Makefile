LIBS_PATH=../libs
CFLAGS=-Wall -Wextra -std=c99

INCS=-I $(LIBS_PATH)/bit_array -I $(LIBS_PATH)/string_buffer \
     -I $(LIBS_PATH)/seq_file -I ../src

LIBS=-L $(LIBS_PATH)/bit_array -L $(LIBS_PATH)/string_buffer -L ../src

LINK=-lalign -lstrbuf -lbitarr -lpthread -lz

all: tests

tests: tests.c ../src/libalign.a
	$(CC) -o $@ $< $(CFLAGS) $(INCS) $(LIBS) $(LINK)

test: tests
	./tests

clean:
	rm -rf tests

.PHONY: all clean test
