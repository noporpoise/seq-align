/*
 needleman_wunsch.c
 project: NeedlemanWunsch
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/needlemanwunsch
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "needleman_wunsch.h"
#include "utility_lib.h"

/*
struct
{
  
} NWAligner;
*/

/* Allocate memory for alignment results */

size_t nw_alloc_mem(const char* seq_a, const char* seq_b,
                 char** alignment_a, char** alignment_b)
{
  size_t length_a = strlen(seq_a);
  size_t length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  size_t longest_alignment = length_a + length_b;
  
  // longest_alignment + 1 to allow for \0
  *alignment_a = (char*) malloc((longest_alignment+1) * sizeof(char));
  *alignment_b = (char*) malloc((longest_alignment+1) * sizeof(char));

  return longest_alignment;
}

// length is = length_a + length_b
char nw_realloc_mem(unsigned int length, char** alignment_a, char** alignment_b)
{
  // longest_alignment + 1 to allow for \0
  char* temp_a = realloc(*alignment_a, (length+1) * sizeof(char));
  char* temp_b = realloc(*alignment_b, (length+1) * sizeof(char));

  if(temp_a != NULL)
  {
    *alignment_a = temp_a;
  }
  
  if(temp_b != NULL)
  {
    *alignment_b = temp_b;
  }

  return (temp_a != NULL && temp_b != NULL);
}

// Find backtrack start when scoring->no_end_gap_penalty is 1
char find_end_max(score_t *score_arr,
                  unsigned int length_a, unsigned int length_b,
                  score_t *curr_score,
                  unsigned int *seq_i, unsigned int *seq_j)
{
  unsigned int i, j;
  score_t temp;
  char updated = 0;

  for(i = 1; i <= length_a; i++)
  {
    temp = ARR_LOOKUP(score_arr, length_a+1, i, length_b);

    if(temp > *curr_score)
    {
      *curr_score = temp;
      *seq_i = i-1;
      *seq_j = length_b-1;
      updated = 1;
    }
  }

  for(j = 1; j <= length_b; j++)
  {
    temp = ARR_LOOKUP(score_arr, length_a+1, length_a, j);

    if(temp > *curr_score)
    {
      *curr_score = temp;
      *seq_i = length_a-1;
      *seq_j = j-1;
      updated = 1;
    }
  }

  return updated;
}

// Alignment between two sequences using needleman-wunsch and a given
// scoring system
int needleman_wunsch(const char* seq_a, const char* seq_b,
                     char* alignment_a, char* alignment_b,
                     const SCORING_SYSTEM* scoring)
{
  #ifdef DEBUG
  scoring_print(scoring);
  #endif

  size_t length_a = strlen(seq_a);
  size_t length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  unsigned int longest_alignment = (unsigned int)length_a +
                                   (unsigned int)length_b;
  
  size_t score_width = length_a+1;
  size_t score_height = length_b+1;
  
  // Check sequences aren't too long to align
  if(score_height > ULONG_MAX / score_width)
  {
    fprintf(stderr, "Error: sequences too long to align (%i * %i > %lu)\n",
            (int)score_width, (int)score_height, ULONG_MAX);
    exit(EXIT_FAILURE);
  }
  
  unsigned long arr_size = (unsigned long)score_width * score_height;
  
  //
  // Check system memory size
  size_t matrix_num_of_bytes = arr_size * sizeof(score_t);
  // 3 matrices required -- overall mem usage:
  size_t overal_num_of_bytes = 3 * matrix_num_of_bytes;

  // Get amount of RAM on this machine
  size_t system_mem_size = getMemorySize();

  if(system_mem_size != 0 && overal_num_of_bytes > system_mem_size)
  {
    fprintf(stderr, "Warning: memory usage exceeds available RAM "
                    "(%lu bytes requested, %lu available)\n",
            overal_num_of_bytes, system_mem_size);
  }

  // 2d array (length_a x length_b)
  // addressing [a][b]

  // Score having just matched
  score_t* match_score = (score_t*) malloc(matrix_num_of_bytes);
  // score having just deleted from seq_a
  score_t* gap_a_score = (score_t*) malloc(matrix_num_of_bytes);
  // score having just inserted into seq_a
  score_t* gap_b_score = (score_t*) malloc(matrix_num_of_bytes);
  
  if(match_score == NULL || gap_a_score == NULL || gap_b_score == NULL)
  {
    fprintf(stderr, "Couldn't allocate enough memory (%lu bytes required)\n",
            matrix_num_of_bytes);
    
    #ifdef DEBUG
    fprintf(stderr, "SeqA length: %i; SeqB length: %i\n",
            (int)length_a, (int)length_b);
    fprintf(stderr, "arr_size: %lu; int size: %li\n", arr_size, sizeof(int));
    #endif

    exit(EXIT_FAILURE);
  }
  
  #ifdef DEBUG
  printf("Malloc'd score matrices!\n");
  #endif

  // Fill in traceback matrix

  if(scoring->no_gaps_in_a)
  {
    unsigned int i;
    for(i = 0; i < arr_size; i++)
    {
      gap_a_score[i] = INT_MIN;
    }
  }

  if(scoring->no_gaps_in_b)
  {
    unsigned int i;
    for(i = 0; i < arr_size; i++)
    {
      gap_b_score[i] = INT_MIN;
    }
  }

  unsigned int i, j;
  
  // [0][0]
  match_score[0] = 0;
  gap_a_score[0] = 0;
  gap_b_score[0] = 0;
  
  // work along first row -> [i][0]
  for(i = 1; i < score_width; i++)
  {
    match_score[i] = INT_MIN;
    
    // Think carefully about which way round these are
    gap_a_score[i] = INT_MIN;
    gap_b_score[i]
      = (score_t)(scoring->no_start_gap_penalty ? 0
                  : scoring->gap_open + (long)i * scoring->gap_extend);
  }

  // work down first column -> [0][j]
  for(j = 1; j < score_height; j++)
  {
    unsigned long index = (unsigned long)j*score_width;
    match_score[index] = INT_MIN;
    
    // Think carefully about which way round these are
    gap_a_score[index]
      = (score_t)(scoring->no_start_gap_penalty ? 0
                  : scoring->gap_open + (long)i * scoring->gap_extend);
    gap_b_score[index] = INT_MIN;
  }
  
  //
  // Update Dynamic Programming arrays
  //

  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  int gap_extend_penalty = scoring->gap_extend;

  for (i = 1; i < score_width; i++)
  {
    for (j = 1; j < score_height; j++)
    {
      // It's an indexing thing...
      unsigned int seq_i = i - 1;
      unsigned int seq_j = j - 1;
      
      // Update match_score[i][j] with position [i-1][j-1]
      // Addressing array must be done with unsigned long
      unsigned long new_index = (unsigned long)j*score_width + i;
      unsigned long old_index;
      
      if(scoring->no_mismatches &&
         !scoring_is_match(scoring, seq_a[seq_i], seq_b[seq_j]))
      {
        match_score[new_index] = INT_MIN;
      }
      else
      {
        old_index = ARR_2D_INDEX(score_width, i-1, j-1);

        // substitution penalty
        int substitution_penalty = scoring_lookup(scoring,
                                                  seq_a[seq_i],
                                                  seq_b[seq_j]);

        // substitution
        match_score[new_index]
          = (score_t)max3((long)match_score[old_index], // continue alignment
                          (long)gap_a_score[old_index], // close gap in seq_a
                          (long)gap_b_score[old_index]) // close gap in seq_b
            + substitution_penalty;
      }                                     

      // Long arithmetic since some INTs are set to INT_MIN and penalty is -ve
      // (adding as ints would cause an integer overflow)

      if(!scoring->no_gaps_in_a)
      {
        // Update gap_a_score[i][j] from position [i][j-1]
        old_index = ARR_2D_INDEX(score_width, i, j-1);

        if(i == score_width-1 && scoring->no_end_gap_penalty)
        {
          gap_a_score[new_index]
            = (score_t)max3((long)match_score[old_index],
                            (long)gap_a_score[old_index],
                            (long)gap_b_score[old_index] + (j == 1 ? 0 : gap_open_penalty));
        }
        else
        {
          gap_a_score[new_index]
            = (score_t)max3((long)match_score[old_index] + gap_open_penalty,
                            (long)gap_a_score[old_index] + gap_extend_penalty,
                            (long)gap_b_score[old_index] + gap_open_penalty);
        }
      }

      if(!scoring->no_gaps_in_b)
      {
        // Update gap_b_score[i][j] from position [i-1][j]
        old_index = ARR_2D_INDEX(score_width, i-1, j);
        
        if(j == score_height-1 && scoring->no_end_gap_penalty)
        {
          gap_b_score[new_index]
            = (score_t)max3((long)match_score[old_index],
                            (long)gap_a_score[old_index] + (i == 1 ? 0 : gap_open_penalty),
                            (long)gap_b_score[old_index]);
        }
        else
        {
          gap_b_score[new_index]
            = (score_t)max3((long)match_score[old_index] + gap_open_penalty,
                            (long)gap_a_score[old_index] + gap_open_penalty,
                            (long)gap_b_score[old_index] + gap_extend_penalty);
        }
      }
    }
  }

  if(scoring->no_gaps_in_a)
  {
    // Allow gaps only at the start/end of A
    unsigned long old_index = ARR_2D_INDEX(score_width, score_width-1, 0);

    for(j = 1; j < score_height; j++)
    {
      unsigned long new_index = ARR_2D_INDEX(score_width, score_width-1, j);

      gap_a_score[new_index]
        = (score_t)max3((long)match_score[old_index] + gap_open_penalty,
                        (long)gap_a_score[old_index] + gap_extend_penalty,
                        INT_MIN);
    
      old_index = new_index;
    }
  }

  if(scoring->no_gaps_in_b)
  {
    // Allow gaps only at the start/end of B
    unsigned long old_index = ARR_2D_INDEX(score_width, 0, score_height-1);

    for(i = 1; i < score_width; i++)
    {
      unsigned long new_index = ARR_2D_INDEX(score_width, i, score_height-1);

      gap_b_score[new_index]
        = (score_t)max3((long)match_score[old_index] + gap_open_penalty,
                        (long)gap_b_score[old_index] + gap_extend_penalty,
                        INT_MIN);
    
      old_index = new_index;
    }
  }

  #ifdef DEBUG
  printf("Filled score matrices - traceback next\n");
  #endif
  
  //
  // Trace back now (score matrices all calculated)
  //
  
  // work backwards re-tracing optimal alignment, then shift sequences into place
  
  // Position of next alignment character in buffer (working backwards)
  unsigned int next_char = longest_alignment-1;

  // Get max score (and therefore current matrix)
  enum Matrix curr_matrix = MATCH;
  score_t curr_score = match_score[arr_size-1];
  
  if(gap_b_score[arr_size-1] >= curr_score)
  {
    curr_matrix = GAP_B;
    curr_score = gap_b_score[arr_size-1];
  }

  if(gap_a_score[arr_size-1] >= curr_score)
  {
    curr_matrix = GAP_A;
    curr_score = gap_a_score[arr_size-1];
  }

  // coords in score matrices
  size_t score_x = length_a;
  size_t score_y = length_b;

#ifdef DEBUG
  alignment_print_matrices(match_score, gap_a_score, gap_b_score,
                           length_a, length_b);
#endif

  // Hold this value to return later
  int max_alignment_score = curr_score;
  
  //unsigned int score_x = (unsigned int)seq_i+1;
  //unsigned int score_y = (unsigned int)seq_j+1;
  unsigned long arr_index;

  // note: longest_alignment = strlen(seq_a) + strlen(seq_b)
  for(; score_x > 0 && score_y > 0; next_char--)
  {
    #ifdef DEBUG
    printf("matrix: %s (%lu,%lu) score: %i\n",
           MATRIX_NAME(curr_matrix), score_x-1, score_y-1, curr_score);
    #endif

    switch(curr_matrix)
    {
      case MATCH:
        alignment_a[next_char] = seq_a[score_x-1];
        alignment_b[next_char] = seq_b[score_y-1];
        break;

      case GAP_A:
        alignment_a[next_char] = '-';
        alignment_b[next_char] = seq_b[score_y-1];
        break;

      case GAP_B:
        alignment_a[next_char] = seq_a[score_x-1];
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
                             &score_x, &score_y, &arr_index,
                             score_width, score_height,
                             // score matrices:
                             match_score, gap_a_score, gap_b_score,
                             seq_a, seq_b, scoring);
    }
  }
  
  // Free memory
  free(match_score);
  free(gap_a_score);
  free(gap_b_score);
  
  // Gap in A
  while(score_y > 0)
  {
    alignment_a[next_char] = '-';
    alignment_b[next_char] = seq_b[score_y-1];
    next_char--;
    score_y--;
  }

  // Gap in B
  while(score_x > 0)
  {
    alignment_a[next_char] = seq_a[score_x-1];
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

  return max_alignment_score;
}
