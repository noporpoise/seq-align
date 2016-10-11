LIBS_PATH=../libs
CFLAGS=-Wall -Wextra -std=c99
CXXFLAGS=-Wall -Wextra

INCS=-I $(LIBS_PATH) -I ../src
LIBS=-L $(LIBS_PATH)/string_buffer -L ../src
LINK=-lalign -lstrbuf -lpthread -lz

all: nw_example sw_example nw_cpp sw_cpp

# C Examples
sw_example: sw_example.c ../src/libalign.a
	$(CC) $< -o $@ $(CFLAGS) $(INCS) $(LIBS) $(LINK)

nw_example: nw_example.c ../src/libalign.a
	$(CC) $< -o $@ $(CFLAGS) $(INCS) $(LIBS) $(LINK)

# C++ Examples
nw_cpp: nw_example.cpp ../src/libalign.a
	$(CXX) $< -o $@ $(CXXFLAGS) $(INCS) $(LIBS) $(LINK)

sw_cpp: sw_example.cpp ../src/libalign.a
	$(CXX) $< -o $@ $(CXXFLAGS) $(INCS) $(LIBS) $(LINK)

clean:
	rm -rf nw_example sw_example nw_cpp sw_cpp *.greg *.dSYM

.PHONY: all clean
