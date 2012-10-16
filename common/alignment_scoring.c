/*
 alignment_scoring.c
 project: AlignmentScoring
 author: Isaac Turner <turner.isaac@gmail.com>
 Used in SmithWaterman and NeedlemanWunsch projects
 url: http://sourceforge.net/projects/needlemanwunsch
 url: http://sourceforge.net/projects/smithwaterman
 Copyright (C) 06-Dec-2011

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

// Turn on debugging output by defining DEBUG
//#define DEBUG

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h> // tolower

#include "alignment_scoring.h"

#define WILDCARDS_INIT_LEN 10

SCORING_SYSTEM* scoring_create(int match, int mismatch,
                               int gap_open, int gap_extend,
                               char no_start_gap_penalty,
                               char no_end_gap_penalty,
                               char no_gaps_in_a, char no_gaps_in_b,
                               char no_mismatches,
                               char case_sensitive)
{
  SCORING_SYSTEM* scoring = (SCORING_SYSTEM*) malloc(sizeof(SCORING_SYSTEM));

  // Gap of length 1 has penalty (gap_open+gap_extend)
  // of length N: (gap_open + gap_extend*N)
  scoring->gap_open = gap_open;
  scoring->gap_extend = gap_extend;

  scoring->no_start_gap_penalty = no_start_gap_penalty;
  scoring->no_end_gap_penalty = no_end_gap_penalty;

  scoring->no_gaps_in_a = no_gaps_in_a;
  scoring->no_gaps_in_b = no_gaps_in_b;
  scoring->no_mismatches = no_mismatches;

  scoring->use_match_mismatch = 1;
  scoring->match = match;
  scoring->mismatch = mismatch;

  scoring->case_sensitive = case_sensitive;

  scoring->num_of_wildcards = 0;
  scoring->wildcards_arrlen = 0;
  scoring->wildcards = NULL;
  scoring->wildscores = NULL;

  scoring->swap_table = NULL;

  return scoring;
}

void scoring_add_wildcard(SCORING_SYSTEM* scoring, char c, int s)
{
  if(scoring->wildcards == NULL)
  {
    // Initialise array
    scoring->wildcards = (char*)malloc(WILDCARDS_INIT_LEN*sizeof(char));
    scoring->wildscores = (int*)malloc(WILDCARDS_INIT_LEN*sizeof(int));
    scoring->wildcards_arrlen = WILDCARDS_INIT_LEN;
    scoring->num_of_wildcards = 0;
  }

  if(scoring->num_of_wildcards == scoring->wildcards_arrlen)
  {
    // Double array capacity
    scoring->wildcards_arrlen *= 2;
    scoring->wildcards = realloc(scoring->wildcards,
                                 scoring->wildcards_arrlen * sizeof(char));
    scoring->wildscores = realloc(scoring->wildscores,
                                  scoring->wildcards_arrlen * sizeof(int));
  }

  if(!scoring->case_sensitive)
    c = tolower(c);

  scoring->wildcards[scoring->num_of_wildcards] = c;
  scoring->wildscores[scoring->num_of_wildcards] = s;
  scoring->num_of_wildcards++;
}

char scoring_check_wildcards(const SCORING_SYSTEM* scoring, char a, char b)
{
  if(!scoring->case_sensitive)
  {
    a = tolower(a);
    b = tolower(b);
  }

  // Check if either characters are wildcards
  unsigned int i;
  for(i = 0; i < scoring->num_of_wildcards; i++)
  {
    if(scoring->wildcards[i] == a || scoring->wildcards[i] == b)
    {
      return 1;
    }
  }

  return 0;
}

char scoring_is_match(const SCORING_SYSTEM* scoring, char a, char b)
{
  if(!scoring->case_sensitive)
  {
    a = tolower(a);
    b = tolower(b);
  }

  return (a == b || scoring_check_wildcards(scoring, a, b));
}

void scoring_add_mutations(SCORING_SYSTEM* scoring,
                           unsigned int num_chars, char* chars, int* scores,
                           char use_match_mismatch)
{
  unsigned int i, j;
  char a, b;
  int score;

  for(i = 0; i < num_chars; i++)
  {
    a = scoring->case_sensitive ? chars[i] : tolower(chars[i]);

    for(j = 0; j < num_chars; j++)
    {
      b = scoring->case_sensitive ? chars[j] : tolower(chars[j]);
      score = ARR_LOOKUP(scores, num_chars, i, j);

      scoring_add_mutation(scoring, a, b, score);
    }
  }

  scoring->use_match_mismatch = use_match_mismatch;
}

void scoring_add_mutation(SCORING_SYSTEM* scoring, char a, char b, int score)
{
  MUTATION_SCORE** hashtable = &(scoring->swap_table);

  MUTATION_SCORE* new_entry = (MUTATION_SCORE*) malloc(sizeof(MUTATION_SCORE));
  new_entry->id = GENERATE_HASH(a, b);
  new_entry->swap_score = score;
  HASH_ADD_INT(*hashtable, id, new_entry);

#ifdef DEBUG
  printf("adding %c -> %c: %i\n", a, b, score);
#endif
}

void scoring_free(SCORING_SYSTEM* scoring)
{
  if(scoring->swap_table != NULL)
  {
    MUTATION_SCORE *curr, *tmp;

    HASH_ITER(hh, scoring->swap_table, curr, tmp) {
      HASH_DEL(scoring->swap_table, curr);
      free(curr);
    }
  }

  free(scoring);
}

void scoring_print(const SCORING_SYSTEM* scoring)
{
  printf("scoring:\n");
  printf("  match: %i; mismatch: %i; (use_match_mismatch: %i)\n",
         scoring->match, scoring->mismatch, scoring->use_match_mismatch);

  printf("  gap_open: %i; gap_extend: %i;\n",
         scoring->gap_open, scoring->gap_extend);

  printf("  no_gaps_in_a: %i; no_gaps_in_b: %i; no_mismatches: %i;\n",
         scoring->no_gaps_in_a, scoring->no_gaps_in_b, scoring->no_mismatches);

  printf("  no_start_gap_penalty: %i; no_end_gap_penalty: %i;\n",
         scoring->no_start_gap_penalty, scoring->no_end_gap_penalty);

  printf("  swap_table: %s\n", (scoring->swap_table == NULL ? "no" : "yes"));
}

int scoring_lookup(const SCORING_SYSTEM* scoring, char a, char b)
//, char* is_match)
{
  if(!scoring->case_sensitive)
  {
    a = tolower(a);
    b = tolower(b);
  }

  //#ifdef DEBUG
  //printf(" scoring_lookup(%c,%c)\n", a, b);
  //#endif

  //*is_match = (a == b);

  // Look up in table
  if(scoring->swap_table != NULL)
  {
    int hash_key = GENERATE_HASH(a,b);
    MUTATION_SCORE* result;

    HASH_FIND_INT(scoring->swap_table, &hash_key, result);

    if(result != NULL)
    {
      return result->swap_score;
    }
  }

  // Check wildcards
  // Wildcards are used in the order they are given
  // e.g. if we specify '--wildcard X 2 --wildcard Y 3' X:Y align with score 2
  unsigned int i;
  for(i = 0; i < scoring->num_of_wildcards; i++)
  {
    if(scoring->wildcards[i] == a || scoring->wildcards[i] == b)
    {
      //*is_match = 1;
      return scoring->wildscores[i];
    }
  }

  // Use match/mismatch
  if(scoring->use_match_mismatch)
  {
    return a == b ? scoring->match : scoring->mismatch;
  }

  // Error
  fprintf(stderr, "Error: Unknown character pair (%c,%c) and "
                   "match/mismatch have not been set\n", a, b);
  exit(EXIT_FAILURE);
}

//
// Some scoring systems
//

// Scoring for protein comparisons of length <35bp
SCORING_SYSTEM* scoring_system_PAM30()
{
  int pam30[576] =
{ 6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2,  0, -1,-13, -8, -2, -3, -3, -3,-17,
 -7,  8, -6,-10, -8, -2, -9, -9, -2, -5, -8,  0, -4, -9, -4, -3, -6, -2,-10, -8, -7, -4, -6,-17,
 -4, -6,  8,  2,-11, -3, -2, -3,  0, -5, -7, -1, -9, -9, -6,  0, -2, -8, -4, -8,  6, -3, -3,-17,
 -3,-10,  2,  8,-14, -2,  2, -3, -4, -7,-12, -4,-11,-15, -8, -4, -5,-15,-11, -8,  6,  1, -5,-17,
 -6, -8,-11,-14, 10,-14,-14, -9, -7, -6,-15,-14,-13,-13, -8, -3, -8,-15, -4, -6,-12,-14, -9,-17,
 -4, -2, -3, -2,-14,  8,  1, -7,  1, -8, -5, -3, -4,-13, -3, -5, -5,-13,-12, -7, -3,  6, -5,-17,
 -2, -9, -2,  2,-14,  1,  8, -4, -5, -5, -9, -4, -7,-14, -5, -4, -6,-17, -8, -6,  1,  6, -5,-17,
 -2, -9, -3, -3, -9, -7, -4,  6, -9,-11,-10, -7, -8, -9, -6, -2, -6,-15,-14, -5, -3, -5, -5,-17,
 -7, -2,  0, -4, -7,  1, -5, -9,  9, -9, -6, -6,-10, -6, -4, -6, -7, -7, -3, -6, -1, -1, -5,-17,
 -5, -5, -5, -7, -6, -8, -5,-11, -9,  8, -1, -6, -1, -2, -8, -7, -2,-14, -6,  2, -6, -6, -5,-17,
 -6, -8, -7,-12,-15, -5, -9,-10, -6, -1,  7, -8,  1, -3, -7, -8, -7, -6, -7, -2, -9, -7, -6,-17,
 -7,  0, -1, -4,-14, -3, -4, -7, -6, -6, -8,  7, -2,-14, -6, -4, -3,-12, -9, -9, -2, -4, -5,-17,
 -5, -4, -9,-11,-13, -4, -7, -8,-10, -1,  1, -2, 11, -4, -8, -5, -4,-13,-11, -1,-10, -5, -5,-17,
 -8, -9, -9,-15,-13,-13,-14, -9, -6, -2, -3,-14, -4,  9,-10, -6, -9, -4,  2, -8,-10,-13, -8,-17,
 -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -7, -6, -8,-10,  8, -2, -4,-14,-13, -6, -7, -4, -5,-17,
  0, -3,  0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5, -6, -2,  6,  0, -5, -7, -6, -1, -5, -3,-17,
 -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -7, -3, -4, -9, -4,  0,  7,-13, -6, -3, -3, -6, -4,-17,
-13, -2, -8,-15,-15,-13,-17,-15, -7,-14, -6,-12,-13, -4,-14, -5,-13, 13, -5,-15,-10,-14,-11,-17,
 -8,-10, -4,-11, -4,-12, -8,-14, -3, -6, -7, -9,-11,  2,-13, -7, -6, -5, 10, -7, -6, -9, -7,-17,
 -2, -8, -8, -8, -6, -7, -6, -5, -6,  2, -2, -9, -1, -8, -6, -6, -3,-15, -7,  7, -8, -6, -5,-17,
 -3, -7,  6,  6,-12, -3,  1, -3, -1, -6, -9, -2,-10,-10, -7, -1, -3,-10, -6, -8,  6,  0, -5,-17,
 -3, -4, -3,  1,-14,  6,  6, -5, -1, -6, -7, -4, -5,-13, -4, -5, -6,-14, -9, -6,  0,  6, -5,-17,
 -3, -6, -3, -5, -9, -5, -5, -5, -5, -5, -6, -5, -5, -8, -5, -3, -4,-11, -7, -5, -5, -5, -5,-17,
-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,  1};

  char *bases = "ARNDCQEGHILKMFPSTWYVBZX*";

  // *->* match: 1
  // *->* mismatch: -17
  // Gap open -9, gap extend -1
  // no_start_gap_penalty, no_end_gap_penalty = 0
  // case_sensitive = 0
  SCORING_SYSTEM* scoring = scoring_create(1, -17, -9, -1, 0, 0, 0, 0, 0, 0);

  // use_match_mismatch=1
  scoring_add_mutations(scoring, 24, bases, pam30, 1);

  return scoring;
}

// Scoring for protein comparisons of length 35-50
SCORING_SYSTEM* scoring_system_PAM70()
{
  int pam70[576] =
{ 5, -4, -2, -1, -4, -2, -1,  0, -4, -2, -4, -4, -3, -6,  0,  1,  1, -9, -5, -1, -1, -1, -2,-11,
 -4,  8, -3, -6, -5,  0, -5, -6,  0, -3, -6,  2, -2, -7, -2, -1, -4,  0, -7, -5, -4, -2, -3,-11,
 -2, -3,  6,  3, -7, -1,  0, -1,  1, -3, -5,  0, -5, -6, -3,  1,  0, -6, -3, -5,  5, -1, -2,-11,
 -1, -6,  3,  6, -9,  0,  3, -1, -1, -5, -8, -2, -7,-10, -4, -1, -2,-10, -7, -5,  5,  2, -3,-11,
 -4, -5, -7, -9,  9, -9, -9, -6, -5, -4,-10, -9, -9, -8, -5, -1, -5,-11, -2, -4, -8, -9, -6,-11,
 -2,  0, -1,  0, -9,  7,  2, -4,  2, -5, -3, -1, -2, -9, -1, -3, -3, -8, -8, -4, -1,  5, -2,-11,
 -1, -5,  0,  3, -9,  2,  6, -2, -2, -4, -6, -2, -4, -9, -3, -2, -3,-11, -6, -4,  2,  5, -3,-11,
  0, -6, -1, -1, -6, -4, -2,  6, -6, -6, -7, -5, -6, -7, -3,  0, -3,-10, -9, -3, -1, -3, -3,-11,
 -4,  0,  1, -1, -5,  2, -2, -6,  8, -6, -4, -3, -6, -4, -2, -3, -4, -5, -1, -4,  0,  1, -3,-11,
 -2, -3, -3, -5, -4, -5, -4, -6, -6,  7,  1, -4,  1,  0, -5, -4, -1, -9, -4,  3, -4, -4, -3,-11,
 -4, -6, -5, -8,-10, -3, -6, -7, -4,  1,  6, -5,  2, -1, -5, -6, -4, -4, -4,  0, -6, -4, -4,-11,
 -4,  2,  0, -2, -9, -1, -2, -5, -3, -4, -5,  6,  0, -9, -4, -2, -1, -7, -7, -6, -1, -2, -3,-11,
 -3, -2, -5, -7, -9, -2, -4, -6, -6,  1,  2,  0, 10, -2, -5, -3, -2, -8, -7,  0, -6, -3, -3,-11,
 -6, -7, -6,-10, -8, -9, -9, -7, -4,  0, -1, -9, -2,  8, -7, -4, -6, -2,  4, -5, -7, -9, -5,-11,
  0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5, -7,  7,  0, -2, -9, -9, -3, -4, -2, -3,-11,
  1, -1,  1, -1, -1, -3, -2,  0, -3, -4, -6, -2, -3, -4,  0,  5,  2, -3, -5, -3,  0, -2, -1,-11,
  1, -4,  0, -2, -5, -3, -3, -3, -4, -1, -4, -1, -2, -6, -2,  2,  6, -8, -4, -1, -1, -3, -2,-11,
 -9,  0, -6,-10,-11, -8,-11,-10, -5, -9, -4, -7, -8, -2, -9, -3, -8, 13, -3,-10, -7,-10, -7,-11,
 -5, -7, -3, -7, -2, -8, -6, -9, -1, -4, -4, -7, -7,  4, -9, -5, -4, -3,  9, -5, -4, -7, -5,-11,
 -1, -5, -5, -5, -4, -4, -4, -3, -4,  3,  0, -6,  0, -5, -3, -3, -1,-10, -5,  6, -5, -4, -2,-11,
 -1, -4,  5,  5, -8, -1,  2, -1,  0, -4, -6, -1, -6, -7, -4,  0, -1, -7, -4, -5,  5,  1, -2,-11,
 -1, -2, -1,  2, -9,  5,  5, -3,  1, -4, -4, -2, -3, -9, -2, -2, -3,-10, -7, -4,  1,  5, -3,-11,
 -2, -3, -2, -3, -6, -2, -3, -3, -3, -3, -4, -3, -3, -5, -3, -1, -2, -7, -5, -2, -2, -3, -3,-11,
-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,  1};
  
  char *bases = "ARNDCQEGHILKMFPSTWYVBZX*";

  // *->* match: 1
  // *->* mismatch: -11
  // Gap open -10, gap extend -1
  // no_start_gap_penalty, no_end_gap_penalty = 0
  // case_sensitive = 0
  SCORING_SYSTEM* scoring = scoring_create(1, -11, -10, -1, 0, 0, 0, 0, 0, 0);

  // use_match_mismatch=1
  scoring_add_mutations(scoring, 24, bases, pam70, 1);

  return scoring;
}

// Scoring for protein comparisons of length 50-85
SCORING_SYSTEM* scoring_system_BLOSUM80()
{
  int blosum80[576]
    = { 7,-3,-3,-3,-1,-2,-2, 0,-3,-3,-3,-1,-2,-4,-1, 2, 0,-5,-4,-1,-3,-2,-1,-8,
       -3, 9,-1,-3,-6, 1,-1,-4, 0,-5,-4, 3,-3,-5,-3,-2,-2,-5,-4,-4,-2, 0,-2,-8,
       -3,-1, 9, 2,-5, 0,-1,-1, 1,-6,-6, 0,-4,-6,-4, 1, 0,-7,-4,-5, 5,-1,-2,-8,
       -3,-3, 2,10,-7,-1, 2,-3,-2,-7,-7,-2,-6,-6,-3,-1,-2,-8,-6,-6, 6, 1,-3,-8,
       -1,-6,-5,-7,13,-5,-7,-6,-7,-2,-3,-6,-3,-4,-6,-2,-2,-5,-5,-2,-6,-7,-4,-8,
       -2, 1, 0,-1,-5, 9, 3,-4, 1,-5,-4, 2,-1,-5,-3,-1,-1,-4,-3,-4,-1, 5,-2,-8,
       -2,-1,-1, 2,-7, 3, 8,-4, 0,-6,-6, 1,-4,-6,-2,-1,-2,-6,-5,-4, 1, 6,-2,-8,
        0,-4,-1,-3,-6,-4,-4, 9,-4,-7,-7,-3,-5,-6,-5,-1,-3,-6,-6,-6,-2,-4,-3,-8,
       -3, 0, 1,-2,-7, 1, 0,-4,12,-6,-5,-1,-4,-2,-4,-2,-3,-4, 3,-5,-1, 0,-2,-8,
       -3,-5,-6,-7,-2,-5,-6,-7,-6, 7, 2,-5, 2,-1,-5,-4,-2,-5,-3, 4,-6,-6,-2,-8,
       -3,-4,-6,-7,-3,-4,-6,-7,-5, 2, 6,-4, 3, 0,-5,-4,-3,-4,-2, 1,-7,-5,-2,-8,
       -1, 3, 0,-2,-6, 2, 1,-3,-1,-5,-4, 8,-3,-5,-2,-1,-1,-6,-4,-4,-1, 1,-2,-8,
       -2,-3,-4,-6,-3,-1,-4,-5,-4, 2, 3,-3, 9, 0,-4,-3,-1,-3,-3, 1,-5,-3,-2,-8,
       -4,-5,-6,-6,-4,-5,-6,-6,-2,-1, 0,-5, 0,10,-6,-4,-4, 0, 4,-2,-6,-6,-3,-8,
       -1,-3,-4,-3,-6,-3,-2,-5,-4,-5,-5,-2,-4,-6,12,-2,-3,-7,-6,-4,-4,-2,-3,-8,
        2,-2, 1,-1,-2,-1,-1,-1,-2,-4,-4,-1,-3,-4,-2, 7, 2,-6,-3,-3, 0,-1,-1,-8,
        0,-2, 0,-2,-2,-1,-2,-3,-3,-2,-3,-1,-1,-4,-3, 2, 8,-5,-3, 0,-1,-2,-1,-8,
       -5,-5,-7,-8,-5,-4,-6,-6,-4,-5,-4,-6,-3, 0,-7,-6,-5,16, 3,-5,-8,-5,-5,-8,
       -4,-4,-4,-6,-5,-3,-5,-6, 3,-3,-2,-4,-3, 4,-6,-3,-3, 3,11,-3,-5,-4,-3,-8,
       -1,-4,-5,-6,-2,-4,-4,-6,-5, 4, 1,-4, 1,-2,-4,-3, 0,-5,-3, 7,-6,-4,-2,-8,
       -3,-2, 5, 6,-6,-1, 1,-2,-1,-6,-7,-1,-5,-6,-4, 0,-1,-8,-5,-6, 6, 0,-3,-8,
       -2, 0,-1, 1,-7, 5, 6,-4, 0,-6,-5, 1,-3,-6,-2,-1,-2,-5,-4,-4, 0, 6,-1,-8,
       -1,-2,-2,-3,-4,-2,-2,-3,-2,-2,-2,-2,-2,-3,-3,-1,-1,-5,-3,-2,-3,-1,-2,-8,
       -8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8,-8, 1};

  char *bases = "ARNDCQEGHILKMFPSTWYVBZX*";

  // *->* match: 1
  // *->* mismatch: -8
  // Gap open -10, gap extend -1
  // no_start_gap_penalty, no_end_gap_penalty = 0
  // case_sensitive = 0
  SCORING_SYSTEM* scoring = scoring_create(1, -8, -10, -1, 0, 0, 0, 0, 0, 0);
  
  // use_match_mismatch=1
  scoring_add_mutations(scoring, 24, bases, blosum80, 1);

  return scoring;
}

// Scoring for protein comparisons of length >85
SCORING_SYSTEM* scoring_system_BLOSUM62()
{
  int blosum62[576]
    = { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4,
       -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4,
       -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4,
       -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4,
        0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,
       -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4,
       -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
        0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4,
       -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4,
       -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4,
       -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4,
       -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4,
       -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4,
       -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4,
       -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4,
        1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4,
        0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4,
       -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4,
       -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4,
        0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4,
       -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4,
       -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
        0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4,
       -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1};
  
  char *bases = "ARNDCQEGHILKMFPSTWYVBZX*";

  // *->* match: 1
  // *->* mismatch: -4
  // Gap open -10, gap extend -1
  // no_start_gap_penalty, no_end_gap_penalty = 0
  // case_sensitive = 0
  SCORING_SYSTEM* scoring = scoring_create(1, -4, -10, -1, 0, 0, 0, 0, 0, 0);
  
  // use_match_mismatch=1
  scoring_add_mutations(scoring, 24, bases, blosum62, 1);

  return scoring;
}

// Scoring system for predicting DNA hybridization
// "Optimization of the BLASTN substitution matrix for prediction of
//   non-specific DNA microarray hybridization" (2009)
// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2831327/
SCORING_SYSTEM* scoring_system_DNA_hybridization()
{
  int sub_matrix[64] = { 2, 2,-4,-4,-4,-4,-4,-4,
                         2, 2,-4,-4,-4,-4,-4,-4,
                        -4,-4, 5, 5,-4,-4,-4,-4,
                        -4,-4, 5, 5,-4,-4,-4,-4,
                        -4,-4,-4,-4, 5, 5,-4,-4,
                        -4,-4,-4,-4, 5, 5,-4,-4,
                        -4,-4,-4,-4,-4,-4, 2, 2,
                        -4,-4,-4,-4,-4,-4, 2, 2};

  char *bases = "AaCcGgTt";

  // match: 0
  // mismatch: 0
  // Gap open -10, gap extend -10
  // no_start_gap_penalty, no_end_gap_penalty = 0
  // case_sensitive = 0
  SCORING_SYSTEM* scoring = scoring_create(0, 0, -10, -10, 0, 0, 0, 0, 0, 0);
  
  // use_match_mismatch=0
  scoring_add_mutations(scoring, 8, bases, sub_matrix, 0);

  return scoring;
}

// Default
SCORING_SYSTEM* scoring_system_default()
{
  int match_default = 1;
  int mismatch_default = -2;
  int gap_open_default = -4;
  int gap_extend_default = -1;

  // no_start_gap_penalty, no_end_gap_penalty = 0
  // case_sensitive = 0
  SCORING_SYSTEM* scoring = scoring_create(match_default, mismatch_default,
                                           gap_open_default, gap_extend_default,
                                           0, 0, 0, 0, 0, 0);

  return scoring;
}
