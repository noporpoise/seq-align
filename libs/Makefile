
all: bit_array string_buffer seq_file
	cd bit_array &&  git pull && make
	cd string_buffer && git pull && make
	cd seq_file && git pull && make

clean:
	cd bit_array && make clean
	cd string_buffer && make clean
	cd seq_file && make clean
	rm -rf BitArray

.PHONY: all clean

bit_array:
	git clone https://github.com/noporpoise/BitArray.git
	mv BitArray bit_array

string_buffer:
	git clone https://github.com/noporpoise/string_buffer.git

seq_file:
	git clone https://github.com/noporpoise/seq_file.git
