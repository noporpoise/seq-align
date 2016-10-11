/*
 smith_waterman.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// Turn on debugging output by defining SEQ_ALIGN_VERBOSE
//#define SEQ_ALIGN_VERBOSE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sort_r/sort_r.h"

#include "smith_waterman.h"
#include "alignment_macros.h"

typedef struct {
  uint32_t *b; size_t l, s; // l is bits, s in uint32_t
} BitSet;

static inline BitSet* bitset_alloc(BitSet *bs, size_t l) {
  bs->s = (l+31)/32;
  bs->l = l;
  bs->b = calloc(sizeof(bs->b[0]), bs->s);
  return bs->b ? bs : NULL;
}

static inline void bitset_dealloc(BitSet *bs) {
  free(bs->b);
  memset(bs, 0, sizeof(*bs));
}

static inline BitSet* bitset_set_length(BitSet *bs, size_t l) {
  size_t ss = (l+31)/32;
  if(ss > bs->s) {
    bs->b = realloc(bs->b, ss);
    memset(bs->b+bs->s, 0, (bs->s-ss)*sizeof(bs->b[0])); // zero new memory
    bs->s = ss;
  }
  bs->l = l;
  return bs->b ? bs : NULL;
}

// For iterating through local alignments
typedef struct
{
  BitSet match_scores_mask;
  size_t *sorted_match_indices, hits_capacity, num_of_hits, next_hit;
} sw_history_t;

// Store alignment here
struct sw_aligner_t
{
  aligner_t aligner;
  sw_history_t history;
};

// Sort indices by their matrix values
// Struct used to pass data to sort_match_indices
typedef struct
{
  score_t *match_scores;
  unsigned int score_width;
} MatrixSort;

// Function passed to sort_r
int sort_match_indices(const void *aa, const void *bb, void *arg)
{
  const size_t *a = aa;
  const size_t *b = bb;
  const MatrixSort *tmp = arg;

  // Recover variables from the struct
  const score_t *match_scores = tmp->match_scores;
  size_t score_width = tmp->score_width;

  long diff = (long)match_scores[*b] - match_scores[*a];

  // Sort by position (from left to right) on seq_a
  if(diff == 0) return (*a % score_width) - (*b % score_width);
  else return diff > 0 ? 1 : -1;
}

static void _init_history(sw_history_t *hist)
{
  bitset_alloc(&hist->match_scores_mask, 256);
  hist->hits_capacity = 256;
  size_t mem = hist->hits_capacity * sizeof(*(hist->sorted_match_indices));
  hist->sorted_match_indices = malloc(mem);
}

static void _ensure_history_capacity(sw_history_t *hist, size_t arr_size)
{
  if(arr_size > hist->hits_capacity) {
    hist->hits_capacity = ROUNDUP2POW(arr_size);
    bitset_set_length(&hist->match_scores_mask, hist->hits_capacity);

    size_t mem = hist->hits_capacity * sizeof(*(hist->sorted_match_indices));
    hist->sorted_match_indices = realloc(hist->sorted_match_indices, mem);
    if(!hist->match_scores_mask.b || !hist->sorted_match_indices) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

sw_aligner_t* smith_waterman_new()
{
  sw_aligner_t *sw = calloc(1, sizeof(sw_aligner_t));
  _init_history(&(sw->history));
  return sw;
}

void smith_waterman_free(sw_aligner_t *sw)
{
  aligner_destroy(&(sw->aligner));
  bitset_dealloc(&sw->history.match_scores_mask);
  free(sw->history.sorted_match_indices);
  free(sw);
}

aligner_t* smith_waterman_get_aligner(sw_aligner_t *sw)
{
  return &sw->aligner;
}

void smith_waterman_align(const char *a, const char *b,
                          const scoring_t *scoring, sw_aligner_t *sw)
{
  smith_waterman_align2(a, b, strlen(a), strlen(b), scoring, sw);
}

void smith_waterman_align2(const char *a, const char *b,
                           size_t len_a, size_t len_b,
                           const scoring_t *scoring, sw_aligner_t *sw)
{
  aligner_t *aligner = &sw->aligner;
  sw_history_t *hist = &sw->history;
  aligner_align(aligner, a, b, len_a, len_b, scoring, 1);

  size_t arr_size = aligner->score_width * aligner->score_height;
  _ensure_history_capacity(hist, arr_size);

  // Set number of hits
  memset(hist->match_scores_mask.b, 0, (hist->match_scores_mask.l+31)/32);
  hist->num_of_hits = hist->next_hit = 0;

  size_t pos;
  for(pos = 0; pos < arr_size; pos++) {
    if(aligner->match_scores[pos] > 0)
      hist->sorted_match_indices[hist->num_of_hits++] = pos;
  }

  // Now sort matched hits
  MatrixSort tmp_struct = {aligner->match_scores, aligner->score_width};
  sort_r(hist->sorted_match_indices, hist->num_of_hits,
         sizeof(size_t), sort_match_indices, &tmp_struct);
}

// Return 1 if alignment was found, 0 otherwise
static char _follow_hit(sw_aligner_t* sw, size_t arr_index,
                        alignment_t* result)
{
  const aligner_t *aligner = &(sw->aligner);
  const sw_history_t *hist = &(sw->history);

  // Follow path through matrix
  size_t score_x = (size_t)ARR_2D_X(arr_index, aligner->score_width);
  size_t score_y = (size_t)ARR_2D_Y(arr_index, aligner->score_width);

  // Local alignments always (start and) end with a match
  enum Matrix curr_matrix = MATCH;
  score_t curr_score = aligner->match_scores[arr_index];

  // Store end arr_index and (x,y) coords for later
  size_t end_arr_index = arr_index;
  size_t end_score_x = score_x;
  size_t end_score_y = score_y;
  score_t end_score = curr_score;

  size_t length;

  for(length = 0; ; length++)
  {
    if(bitset32_get(hist->match_scores_mask.b, arr_index)) return 0;
    bitset32_set(hist->match_scores_mask.b, arr_index);

    if(curr_score == 0) break;

    //printf(" %i (%i, %i)\n", length, score_x, score_y);

    // Find out where to go next
    alignment_reverse_move(&curr_matrix, &curr_score,
                           &score_x, &score_y, &arr_index, aligner);
  }

  // We got a result!
  // Allocate memory for the result
  result->length = length;

  alignment_ensure_capacity(result, length);

  // Jump back to the end of the alignment
  arr_index = end_arr_index;
  score_x = end_score_x;
  score_y = end_score_y;
  curr_matrix = MATCH;
  curr_score = end_score;

  // Now follow backwards again to create alignment!
  unsigned int i;

  for(i = length-1; curr_score > 0; i--)
  {
    switch(curr_matrix)
    {
      case MATCH:
        result->result_a[i] = aligner->seq_a[score_x-1];
        result->result_b[i] = aligner->seq_b[score_y-1];
        break;

      case GAP_A:
        result->result_a[i] = '-';
        result->result_b[i] = aligner->seq_b[score_y-1];
        break;

      case GAP_B:
        result->result_a[i] = aligner->seq_a[score_x-1];
        result->result_b[i] = '-';
        break;

      default:
      fprintf(stderr, "Program error: invalid matrix in _follow_hit()\n");
      fprintf(stderr, "Please submit a bug report to: turner.isaac@gmail.com\n");
      exit(EXIT_FAILURE);
    }

    alignment_reverse_move(&curr_matrix, &curr_score,
                           &score_x, &score_y, &arr_index, aligner);
  }

  result->result_a[length] = '\0';
  result->result_b[length] = '\0';

  result->score = end_score;

  result->pos_a = score_x;
  result->pos_b = score_y;

  result->len_a = end_score_x - score_x;
  result->len_b = end_score_y - score_y;

  return 1;
}

int smith_waterman_fetch(sw_aligner_t *sw, alignment_t *result)
{
  sw_history_t *hist = &(sw->history);

  while(hist->next_hit < hist->num_of_hits)
  {
    size_t arr_index = hist->sorted_match_indices[hist->next_hit++];
    // printf("hit %lu/%lu\n", hist->next_hit, hist->num_of_hits);

    if(!bitset32_get(hist->match_scores_mask.b, arr_index) &&
       _follow_hit(sw, arr_index, result))
    {
      return 1;
    }
  }

  return 0;
}
