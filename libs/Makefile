
all: bit_array string_buffer htslib seq_file
	cd htslib; git pull; make
	cd bit_array; git pull; make
	cd string_buffer; git pull; make
	cd seq_file; git pull; make HTSLIB=../htslib

clean:
	cd htslib; make clean
	cd bit_array; make clean
	cd string_buffer; make clean
	cd seq_file; make clean
	rm -rf BitArray

.PHONY: all clean

htslib:
	git clone git://github.com/samtools/htslib.git

bit_array:
	git clone git@github.com:noporpoise/BitArray.git
	mv BitArray bit_array

string_buffer:
	git clone git@github.com:noporpoise/string_buffer.git

seq_file:
	git clone git@github.com:noporpoise/seq_file.git
