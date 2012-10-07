/*
 alignment_scoring.h
 project: AlignmentScoring
 author: Isaac Turner <turner.isaac@gmail.com>
 Used in SmithWaterman and NeedlemanWunsch projects
 url: http://sourceforge.net/projects/needlemanwunsch
 url: http://sourceforge.net/projects/smithwaterman
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

#ifndef ALIGNMENT_SCORING_HEADER_SEEN
#define ALIGNMENT_SCORING_HEADER_SEEN

#define ARR_2D_INDEX(width,i,j) (((unsigned long)(j)*(width)) + (i))
#define ARR_LOOKUP(arr,width,i,j) arr[ARR_2D_INDEX((width),(i),(j))]
#define ARR_2D_X(arr_index, arr_width) ((arr_index) % (arr_width))
#define ARR_2D_Y(arr_index, arr_width) ((arr_index) / (arr_width))

#define GENERATE_HASH(a,b) (int)(((int)(a) << 8) | (int)(b))

#include "uthash.h"

typedef struct MUTATION_SCORE MUTATION_SCORE;
typedef struct SCORING_SYSTEM SCORING_SYSTEM;

struct MUTATION_SCORE
{
  int id; // hash key
  int swap_score;
  UT_hash_handle hh; // makes this structure hashable
};

struct SCORING_SYSTEM
{
  int gap_open, gap_extend;

  // Needleman Wunsch only
  // Turn these on to turn off penalties for gaps at the start/end of alignment
  char no_start_gap_penalty;
  char no_end_gap_penalty;

  // Turn at most one of these on at a time to prevent gaps/mismatches
  char no_gaps;
  char no_mismatches;

  // If swap_table != NULL, but char->char pair swap is not in the hashtable,
  // should we use match/mismatch values?
  char use_match_mismatch;
  int match, mismatch;

  char case_sensitive;

  // Array of characters that match to everything with the same penalty (i.e. 'N's)
  unsigned int num_of_wildcards, wildcards_arrlen;
  char *wildcards;
  int *wildscores;

  MUTATION_SCORE* swap_table;
};

// Constructor
SCORING_SYSTEM* scoring_create(int match, int mismatch,
                               int gap_open, int gap_extend,
                               char no_start_gap_penalty,
                               char no_end_gap_penalty,
                               char no_gaps, char no_mismatches,
                               char case_sensitive);

void scoring_add_wildcard(SCORING_SYSTEM* scoring, char c, int s);

char scoring_check_wildcards(const SCORING_SYSTEM* scoring, char a, char b);

void scoring_add_mutations(SCORING_SYSTEM* scoring,
                           unsigned int num_chars, char* chars, int* scores,
                           char use_match_mismatch);

void scoring_add_mutation(SCORING_SYSTEM* scoring, char a, char b, int score);
void scoring_free(SCORING_SYSTEM* scoring);

void scoring_print(const SCORING_SYSTEM* scoring);

int scoring_lookup(const SCORING_SYSTEM* scoring, char a, char b);

// Some scoring systems
SCORING_SYSTEM* scoring_system_PAM30();
SCORING_SYSTEM* scoring_system_PAM70();
SCORING_SYSTEM* scoring_system_BLOSUM80();
SCORING_SYSTEM* scoring_system_BLOSUM62();
SCORING_SYSTEM* scoring_system_DNA_hybridization();
SCORING_SYSTEM* scoring_system_default(); // DNA/RNA default

#endif
