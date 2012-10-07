# set this
ifndef SWPATH
	SWPATH:=../
endif

LIBS_PATH := ../../libs
ALIGN_PATH := $(SWPATH)/$(LIBS_PATH)/alignment_scoring/
BIT_ARRAY_PATH := $(SWPATH)/$(LIBS_PATH)/bit_array/

FLAGS := -DSCORE_TYPE='unsigned int' -I $(SWPATH) -I $(ALIGN_PATH) -I $(BIT_ARRAY_PATH)

# Check mac/linux
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	FLAGS := $(FLAGS) -fnested-functions
endif

all:
	gcc -o sw_example -Wall $(FLAGS) $(SWPATH)/smith_waterman.c $(ALIGN_PATH)/alignment.c $(ALIGN_PATH)/alignment_scoring.c  $(BIT_ARRAY_PATH)/bit_array.c sw_example.c

clean:
	if test -e sw_example; then rm sw_example; fi
	if test -e sw_example.dSYM; then rm -r sw_example.dSYM; fi
	if test -e sw_example.greg; then rm sw_example.greg; fi
