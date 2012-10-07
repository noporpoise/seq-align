# set this
ifndef NWPATH
	NWPATH:=../
endif

LIBS_PATH=../../libs
ALIGN_PATH:=$(NWPATH)/$(LIBS_PATH)/alignment_scoring/

FLAGS := -DSCORE_TYPE='int' -I $(NWPATH) -I $(ALIGN_PATH)

all:
	gcc -o nw_example -Wall $(FLAGS) $(NWPATH)/needleman_wunsch.c $(ALIGN_PATH)/alignment.c $(ALIGN_PATH)/alignment_scoring.c  nw_example.c

clean:
	if test -e nw_example; then rm nw_example; fi
	if test -e nw_example.dSYM; then rm -r nw_example.dSYM; fi
	if test -e nw_example.greg; then rm nw_example.greg; fi
