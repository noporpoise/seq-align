/*
 sw_example.c
 project: SmithWaterman
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/smithwaterman/
 Copyright (C) 30-Jan-2012

 see README
 
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

#include "smith_waterman.h"

void align(char* seq_a, char* seq_b, int min_score)
{
  // Decide on scoring
  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  // Compare character case-sensitively (usually set to 0 for DNA etc)
  char case_sensitive = 0;

  // Create scoring system
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

  // Do alignment
  SW_COMPUTATION* smithwaterman = smith_waterman_align(seq_a, seq_b, scoring);

  // Allocate memory for storing result
  SW_LOCAL_ALIGNMENT* alignment
    = (SW_LOCAL_ALIGNMENT*) malloc(sizeof(SW_LOCAL_ALIGNMENT));
  
  // Loop through results
  while(smith_waterman_get_hit(smithwaterman, alignment) &&
        alignment->score >= min_score)
  {
    printf("seqA [%u]: %s\n", alignment->pos_a, alignment->result_a);
    printf("seqB [%u]: %s\n", alignment->pos_b, alignment->result_b);
    printf("score: %i\n\n", alignment->score);
  }

  // Free result
  free(alignment);

  // Free memory used to store scoring preferences
  scoring_free(scoring);
}

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    printf("usage: ./sw_example <seq1> <seq2>\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    align(argv[1], argv[2], 2);
  }

  return 0;
}
