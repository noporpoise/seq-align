/*
 alignment.c
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

#include "alignment.h"

char* align_col_mismatch = "\033[92m"; // Mismatch (GREEN)
char* align_col_indel = "\033[91m"; // Insertion / deletion (RED)
// Pink used by SmithWaterman local alignment for printing surrounding bases
char* align_col_context = "\033[95m";
char* align_col_stop = "\033[0m";

/*
struct Alignment
{
  // Store local alignment result here
  char *result_a, *result_b;
  unsigned int capacity, length;
  unsigned int pos_a, pos_b; // position of first base (0-based)
  unsigned int len_a, len_b; // number of bases in alignment
  score_t score;
};
*/

long max2(long a, long b)
{
  return (a >= b ? a : b);
}

long max3(long a, long b, long c)
{
  long result = a;

  if(b > result) {
    result = b;
  }
  if(c > result) {
    result = c;
  }

  return result;
}

long max4(long a, long b, long c, long d)
{
  long result = a;

  if(b > result) {
    result = b;
  }
  if(c > result) {
    result = c;
  }
  if(d > result) {
    result = d;
  }
  
  return result;
}


void alignment_print_matrices(const score_t* match_score,
                              const score_t* gap_a_score,
                              const score_t* gap_b_score,
                              int length_a, int length_b)
{
  int score_width = length_a+1;
  int i, j;

  printf("match_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%3i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", (int)ARR_LOOKUP(match_score, score_width, i, j));
    }
    printf("\n");
  }
  printf("gap_a_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%3i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", (int)ARR_LOOKUP(gap_a_score, score_width, i, j));
    }
    printf("\n");
  }
  printf("gap_b_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%3i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", (int)ARR_LOOKUP(gap_b_score, score_width, i, j));
    }
    printf("\n");
  }
}

void alignment_colour_print_against(const char *alignment_a,
                                    const char *alignment_b,
                                    char case_sensitive)
{
  int i;
  char red = 0, green = 0;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_b[i] == '-')
    {
      if(!red)
      {
        printf("%s", align_col_indel);
        red = 1;
      }
    }
    else if(red)
    {
      red = 0;
      printf("%s", align_col_stop);
    }

    if(((case_sensitive && alignment_a[i] != alignment_b[i]) ||
        (!case_sensitive && tolower(alignment_a[i]) != tolower(alignment_b[i]))) &&
       alignment_a[i] != '-' && alignment_b[i] != '-')
    {
      if(!green)
      {
        printf("%s", align_col_mismatch);
        green = 1;
      }
    }
    else if(green)
    {
      green = 0;
      printf("%s", align_col_stop);
    }

    printf("%c", alignment_a[i]);
  }

  if(green || red)
  {
    // Stop all colours
    printf("%s", align_col_stop);
  }
}

// Order of alignment_a / alignment_b is not important
void alignment_print_spacer(const char* alignment_a, const char* alignment_b,
                            const SCORING_SYSTEM* scoring)
{
  int i;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-' || alignment_b[i] == '-')
    {
      printf(" ");
    }
    else if(alignment_a[i] == alignment_b[i] ||
            (!scoring->case_sensitive &&
             tolower(alignment_a[i]) == tolower(alignment_b[i])))
    {
      printf("|");
    }
    else
    {
      printf("*");
    }
  }
}

// Backtrack through scoring matrices
void alignment_reverse_move(enum Matrix *curr_matrix, score_t* curr_score,
                            size_t *score_x, size_t *score_y,
                            unsigned long *arr_index,
                            size_t score_width, size_t score_height,
                            const score_t *match_score,
                            const score_t *gap_a_score,
                            const score_t *gap_b_score,
                            const char* seq_a, const char* seq_b,
                            const SCORING_SYSTEM* scoring)
{
  int prev_match_penalty, prev_gap_a_penalty, prev_gap_b_penalty;

  size_t seq_x = (*score_x)-1;
  size_t seq_y = (*score_y)-1;

  int match_penalty = scoring_lookup(scoring, seq_a[seq_x], seq_b[seq_y]);

  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  int gap_extend_penalty = scoring->gap_extend;

  if(scoring->no_end_gap_penalty &&
     (*score_x == score_width-1 || *score_y == score_height-1))
  {
    gap_open_penalty = 0;
    gap_extend_penalty = 0;
  }

  switch(*curr_matrix)
  {
    case MATCH:
      prev_match_penalty = match_penalty;
      prev_gap_a_penalty = match_penalty;
      prev_gap_b_penalty = match_penalty;
      (*score_x)--;
      (*score_y)--;
      break;

    case GAP_A:
      prev_match_penalty = gap_open_penalty;
      prev_gap_a_penalty = gap_extend_penalty;
      prev_gap_b_penalty = gap_open_penalty;
      (*score_y)--;
      break;

    case GAP_B:
      prev_match_penalty = gap_open_penalty;
      prev_gap_a_penalty = gap_open_penalty;
      prev_gap_b_penalty = gap_extend_penalty;
      (*score_x)--;
      break;

    default:
      fprintf(stderr, "Program error: invalid matrix in get_reverse_move()\n");
      fprintf(stderr, "Please submit a bug report to: turner.isaac@gmail.com\n");
      exit(EXIT_FAILURE);
  }

  *arr_index = ARR_2D_INDEX(score_width, *score_x, *score_y);

  if((!scoring->no_gaps_in_a || *score_x == 0 || *score_x == score_width-1) &&
     (long)gap_a_score[*arr_index] + prev_gap_a_penalty == *curr_score)
  {
    *curr_matrix = GAP_A;
    *curr_score = gap_a_score[*arr_index];
    return;
  }
  else if((!scoring->no_gaps_in_b || *score_y == 0 || *score_y == score_height-1) &&
          (long)gap_b_score[*arr_index] + prev_gap_b_penalty == *curr_score)
  {
    *curr_matrix = GAP_B;
    *curr_score = gap_b_score[*arr_index];
  }
  else if(!scoring->no_mismatches &&
          (long)match_score[*arr_index] + prev_match_penalty == *curr_score)
  {
    *curr_matrix = MATCH;
    *curr_score = match_score[*arr_index];
  }
  else
  {
    fprintf(stderr, "Program error: traceback fail (get_reverse_move)\n");
    fprintf(stderr, "This may be due to an integer overflow if your "
                    "sequences are long or if scores are large.  \n");
    fprintf(stderr, "If this is the case using smaller scores or "
                    "shorter sequences may work around this problem.  \n");
    fprintf(stderr, " If you think this is a bug, please report it to: "
                    "turner.isaac@gmail.com\n");
    exit(EXIT_FAILURE);
  }
}
