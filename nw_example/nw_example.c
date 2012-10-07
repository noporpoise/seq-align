/*
 nw_example.c
 project: NeedlemanWunsch
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 06-Dec-2011
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "needleman_wunsch.h"

void align(char* seq_a, char* seq_b)
{
  // Variables to store alignment result
  char *alignment_a, *alignment_b;

  // malloc the above variables
  // (seq1 and seq2 are used to figure out how much memory may be needed)
  nw_alloc_mem(seq_a, seq_b, &alignment_a, &alignment_b);

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

  // Compare character case-sensitively (usually set to 0 for DNA etc)
  char case_sensitive = 0;

  SCORING_SYSTEM* scoring = scoring_create(match, mismatch,
                                           gap_open, gap_extend,
                                           no_start_gap_penalty,
                                           no_end_gap_penalty,
                                           case_sensitive);

  // Add some special cases
  // x -> y means x in seq1 changing to y in seq2
  scoring_add_mutation(scoring, 'a', 'c', -2); // a -> c give substitution score -2
  scoring_add_mutation(scoring, 'c', 'a', -1); // c -> a give substitution score -1

  // We could also prohibit the aligning of characters not given as special cases
  // scoring->use_match_mismatch = 0;

  int score = needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

  printf("seqA: %s\n", alignment_a);
  printf("seqB: %s\n", alignment_b);
  printf("alignment score: %i\n", score);

  // Free memory used to store scoring preferences
  scoring_free(scoring);

  free(alignment_a);
  free(alignment_b);
}

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    printf("usage: ./nw_example <seq1> <seq2>\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    align(argv[1], argv[2]);
  }

  return 0;
}
