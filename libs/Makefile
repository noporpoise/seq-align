
TGTS=bit_array string_buffer seq_file sort_r

all: $(TGTS)

#
# Update (make)
#

string_buffer:
	cd string_buffer && $(MAKE)

seq_file:
	cd seq_file && $(MAKE)

sort_r:
	cd sort_r && $(MAKE)

clean:
	[ -d bit_array ] && cd bit_array && $(MAKE) clean
	[ -d string_buffer ] && cd string_buffer && $(MAKE) clean
	[ -d seq_file ] && cd seq_file && $(MAKE) clean
	[ -d sort_r ] && cd sort_r && $(MAKE) clean

.PHONY: all clean $(TGTS)
