/*
 needleman_wunsch.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "needleman_wunsch.h"

nw_aligner_t* needleman_wunsch_new()
{
  nw_aligner_t *nw = calloc(1, sizeof(nw_aligner_t));
  return nw;
}

void needleman_wunsch_free(nw_aligner_t *nw)
{
  aligner_destroy(nw);
  free(nw);
}

void needleman_wunsch_align(const char *a, const char *b,
                            const scoring_t *scoring,
                            nw_aligner_t *nw, alignment_t *result)
{
  needleman_wunsch_align2(a, b, strlen(a), strlen(b), scoring, nw, result);
}

void needleman_wunsch_align2(const char *a, const char *b,
                             size_t len_a, size_t len_b,
                             const scoring_t *scoring,
                             nw_aligner_t *nw, alignment_t *result)
{
  aligner_align(nw, a, b, len_a, len_b, scoring, 0);

  // work backwards re-tracing optimal alignment, then shift sequences into place

  // note: longest_alignment = strlen(seq_a) + strlen(seq_b)
  size_t longest_alignment = nw->score_width-1 + nw->score_height-1;
  alignment_ensure_capacity(result, longest_alignment);

  // Position of next alignment character in buffer (working backwards)
  size_t next_char = longest_alignment-1;

  size_t arr_size = nw->score_width * nw->score_height;

  // Get max score (and therefore current matrix)
  enum Matrix curr_matrix = MATCH;
  score_t curr_score = nw->match_scores[arr_size-1];

  if(nw->gap_b_scores[arr_size-1] >= curr_score)
  {
    curr_matrix = GAP_B;
    curr_score = nw->gap_b_scores[arr_size-1];
  }

  if(nw->gap_a_scores[arr_size-1] >= curr_score)
  {
    curr_matrix = GAP_A;
    curr_score = nw->gap_a_scores[arr_size-1];
  }

  #ifdef SEQ_ALIGN_VERBOSE
    alignment_print_matrices(nw);
  #endif

  result->score = curr_score;
  char *alignment_a = result->result_a, *alignment_b = result->result_b;

  // coords in score matrices
  size_t score_x = nw->score_width-1, score_y = nw->score_height-1;
  size_t arr_index = arr_size - 1;

  for(; score_x > 0 && score_y > 0; next_char--)
  {
    #ifdef SEQ_ALIGN_VERBOSE
    printf("matrix: %s (%lu,%lu) score: %i\n",
           MATRIX_NAME(curr_matrix), score_x-1, score_y-1, curr_score);
    #endif

    switch(curr_matrix)
    {
      case MATCH:
        alignment_a[next_char] = nw->seq_a[score_x-1];
        alignment_b[next_char] = nw->seq_b[score_y-1];
        break;

      case GAP_A:
        alignment_a[next_char] = '-';
        alignment_b[next_char] = nw->seq_b[score_y-1];
        break;

      case GAP_B:
        alignment_a[next_char] = nw->seq_a[score_x-1];
        alignment_b[next_char] = '-';
        break;

      default:
        fprintf(stderr, "Program error: invalid matrix number\n");
        fprintf(stderr, "Please submit a bug report to: turner.isaac@gmail.com\n");
        exit(EXIT_FAILURE);
    }

    if(score_x > 0 && score_y > 0)
    {
      alignment_reverse_move(&curr_matrix, &curr_score,
                             &score_x, &score_y, &arr_index, nw);
    }
  }

  // Gap in A
  while(score_y > 0)
  {
    alignment_a[next_char] = '-';
    alignment_b[next_char] = nw->seq_b[score_y-1];
    next_char--;
    score_y--;
  }

  // Gap in B
  while(score_x > 0)
  {
    alignment_a[next_char] = nw->seq_a[score_x-1];
    alignment_b[next_char] = '-';
    next_char--;
    score_x--;
  }

  // Shift alignment strings back into 0th position in char arrays
  int first_char = next_char+1;
  int alignment_len = longest_alignment - first_char;

  // Use memmove
  memmove(alignment_a, alignment_a+first_char, alignment_len);
  memmove(alignment_b, alignment_b+first_char, alignment_len);

  alignment_a[alignment_len] = '\0';
  alignment_b[alignment_len] = '\0';

  result->length = alignment_len;
}
