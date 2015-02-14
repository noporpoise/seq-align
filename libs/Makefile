
TGTS=bit_array string_buffer seq_file sort_r

all: $(TGTS)

#
# Download (git clone)
#

bit_array/Makefile:
	git clone https://github.com/noporpoise/BitArray.git bit_array

string_buffer/Makefile:
	git clone https://github.com/noporpoise/string_buffer.git

seq_file/Makefile:
	git clone https://github.com/noporpoise/seq_file.git

sort_r/Makefile:
	git clone https://github.com/noporpoise/sort_r.git

#
# Update (git pull + make)
#

bit_array: bit_array/Makefile
	cd bit_array && git pull && make

string_buffer: string_buffer/Makefile
	cd string_buffer && git pull && make

seq_file: seq_file/Makefile
	cd seq_file && git pull && make

sort_r: sort_r/Makefile
	cd sort_r && git pull && make

clean:
	[ -d bit_array ] && cd bit_array && make clean
	[ -d string_buffer ] && cd string_buffer && make clean
	[ -d seq_file ] && cd seq_file && make clean
	[ -d sort_r ] && cd sort_r && make clean

.PHONY: all clean $(TGTS)
