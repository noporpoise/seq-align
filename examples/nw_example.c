/*
 examples/nw_example.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: May 2013
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "needleman_wunsch.h"

void align(char* seq_a, char* seq_b)
{
  // Variables to store alignment result
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  // Decide on scoring
  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  
  // Don't penalise gaps at the start
  // ACGATTT
  // ----TTT would score +3 (when match=+1)
  char no_start_gap_penalty = 1;
  
  // ..or gaps at the end e.g.
  // ACGATTT
  // ACGA--- would score +4 (when match=+1)
  char no_end_gap_penalty = 1;

  char no_gaps_in_a = 0, no_gaps_in_b = 0;
  char no_mismatches = 0;

  // Compare character case-sensitively (usually set to 0 for DNA etc)
  char case_sensitive = 0;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);

  // Add some special cases
  // x -> y means x in seq1 changing to y in seq2
  scoring_add_mutation(&scoring, 'a', 'c', -2); // a -> c give substitution score -2
  scoring_add_mutation(&scoring, 'c', 'a', -1); // c -> a give substitution score -1

  // We could also prohibit the aligning of characters not given as special cases
  // scoring.use_match_mismatch = 0;

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);

  printf("seqA: %s\n", result->result_a);
  printf("seqB: %s\n", result->result_b);
  printf("alignment score: %i\n", result->score);

  // Free memory for storing alignment results
  needleman_wunsch_free(nw);
  alignment_free(result);
}

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    printf("usage: ./nw_example <seq1> <seq2>\n");
    exit(EXIT_FAILURE);
  }

  align(argv[1], argv[2]);
  exit(EXIT_SUCCESS);
}
