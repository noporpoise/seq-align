/*
 needleman_wunsch.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */


#ifndef NEEDLEMAN_WUNSCH_HEADER_SEEN
#define NEEDLEMAN_WUNSCH_HEADER_SEEN

#include "alignment.h"

typedef aligner_t nw_aligner_t;

nw_aligner_t* needleman_wunsch_new();
void needleman_wunsch_free(nw_aligner_t *nw);

void needleman_wunsch_align(const char *a, const char *b,
                            const scoring_t *scoring,
                            nw_aligner_t *nw, alignment_t *result);


#endif
