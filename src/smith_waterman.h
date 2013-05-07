/*
 smith_waterman.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#ifndef SMITH_WATERMAN_HEADER_SEEN
#define SMITH_WATERMAN_HEADER_SEEN

#include "bit_array.h"
#include "alignment.h"

// Store alignment here
typedef struct
{
  aligner_t aligner;

  // For iterating through local alignments
  BIT_ARRAY* match_scores_mask;
  size_t *sorted_match_indices, hits_capacity, num_of_hits, next_hit;
} sw_aligner_t;

sw_aligner_t *smith_waterman_new();
void smith_waterman_free(sw_aligner_t *sw_aligner);

/*
 Do not alter seq_a, seq_b or scoring whilst calling this method
 or between calls to smith_waterman_get_hit
*/
void smith_waterman_align(const char *seq_a, const char *seq_b,
                          const scoring_t *scoring, sw_aligner_t *sw);

// An alignment to read from, and a pointer to memory to store the result
// returns 1 if an alignment was read, 0 otherwise
int smith_waterman_fetch(sw_aligner_t *sw, alignment_t *result);

#endif
