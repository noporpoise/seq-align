/*
 smith_waterman.c
 project: SmithWaterman
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/smithwaterman/
 Copyright (C) 18-Dec-2011

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
#include <string.h>
#include <limits.h>

#include "smith_waterman.h"
#include "utility_lib.h"
#include "mem_size.h"

// Store alignment here
struct SW_COMPUTATION
{
  // keep pointer to scoring system
  SCORING_SYSTEM* scoring;

  // Pointers to aligned sequence
  const char *seq_a, *seq_b;

  // We allocate everything below here
  score_t* match_score;
  score_t* gap_a_score;
  score_t* gap_b_score;

  // For iterating through local alignments
  BIT_ARRAY* match_scores_mask;
  unsigned long *sorted_match_indices;
  unsigned long num_of_hits;

  unsigned long next_hit;

  size_t score_width, score_height;
};


size_t smith_waterman_seq_a_strlen(SW_COMPUTATION *sw)
{
  return sw->score_width-1;
}

size_t smith_waterman_seq_b_strlen(SW_COMPUTATION *sw)
{
  return sw->score_width-1;
}

#ifdef DEBUG

void print_computation(SW_COMPUTATION* comp)
{
  printf("Smith-Waterman COMPUTATION:\n");

  char* str = bit_array_to_string(comp->match_scores_mask);
  printf(" bit mask: %s\n", str);
  free(str);

  unsigned int i;
  printf(" hits: %lu", comp->sorted_match_indices[0]);
  for(i = 1; i < comp->num_of_hits; i++)
  {
    printf(",%lu", comp->sorted_match_indices[i]);
  }
  printf("\n");
}

#endif


// Sort indices by their matrix values
// Struct used to pass data to sort_match_indices
typedef struct
{
  score_t *match_score;
  unsigned int score_width;
} MatrixSort;

// Function passed to sort_r
int sort_match_indices(const void *aa, const void *bb, void *arg)
{
  const int *a = aa;
  const int *b = bb;
  const MatrixSort *tmp = arg;

  // Recover variables from the struct
  const score_t *match_score = tmp->match_score;
  unsigned int score_width = tmp->score_width;

  long diff = (long)match_score[*b] - match_score[*a];

  if(diff == 0)
  {
    // Sort by position (from left to right) on seq_a
    return (*a % score_width) - (*b % score_width);
  }
  else
  {
    return diff > 0 ? 1 : -1;
  }
}

/*
 Do not alter seq_a, seq_b or scoring whilst calling this method,
 or between calls to smith_waterman_get_hit
*/
SW_COMPUTATION* smith_waterman_align(const char* seq_a, const char* seq_b,
                                     SCORING_SYSTEM* scoring)
{
  #ifdef DEBUG
  scoring_print(scoring);
  #endif

  size_t length_a = strlen(seq_a);
  size_t length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  //int longest_alignment = length_a + length_b;
  
  unsigned int score_width = (unsigned int)length_a+1;
  unsigned int score_height = (unsigned int)length_b+1;
  
  // Check sequences aren't too long to align
  if(score_height > ULONG_MAX / score_width)
  {
    fprintf(stderr, "Error: sequences too long to align (%i * %i > %lu)\n",
            score_width, score_height, ULONG_MAX);
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
  score_t* match_score = (score_t*) malloc(arr_size * sizeof(score_t));
  // score having just deleted from seq_a
  score_t* gap_a_score = (score_t*) malloc(arr_size * sizeof(score_t));
  // score having just inserted into seq_a
  score_t* gap_b_score = (score_t*) malloc(arr_size * sizeof(score_t));
  
  if(match_score == NULL || gap_a_score == NULL || gap_b_score == NULL)
  {
    unsigned long num_of_bytes = arr_size * sizeof(int);
    fprintf(stderr, "Couldn't allocate enough memory (%lu bytes required)\n",
            num_of_bytes);
    
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
  
  unsigned int i, j;
  
  // [0][0]
  match_score[0] = 0;
  gap_a_score[0] = 0;
  gap_b_score[0] = 0;
  
  // work along first row -> [i][0]
  for(i = 1; i < score_width; i++)
  {
    match_score[i] = 0;
    gap_a_score[i] = 0;
    gap_b_score[i] = 0;
  }
  
  // work down first column -> [0][j]
  for(j = 1; j < score_height; j++)
  {
    unsigned long index = (unsigned long)j*score_width;
    match_score[index] = 0;
    gap_a_score[index] = 0;
    gap_b_score[index] = 0;
  }

  if(scoring->no_gaps)
  {
    memset(gap_a_score, 0, matrix_num_of_bytes);
    memset(gap_b_score, 0, matrix_num_of_bytes);
  }

  //memset(match_score, 0, matrix_num_of_bytes);

  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;

  for(i = 1; i < score_width; i++)
  {
    for(j = 1; j < score_height; j++)
    {
      // It's an indexing thing...
      int seq_i = i - 1;
      int seq_j = j - 1;

      // Addressing array must be done with unsigned long
      unsigned long new_index = ARR_2D_INDEX(score_width, i, j);
      unsigned long old_index;

      // 1) Update match_score[i][j] from position [i-1][j-1]

      if(scoring->no_mismatches &&
         !scoring_is_match(scoring, seq_a[seq_i], seq_b[seq_j]))
      {
        match_score[new_index] = 0;
      }
      else
      {
        old_index = ARR_2D_INDEX(score_width, i-1, j-1);

        //substitution penalty
        int substitution_penalty = scoring_lookup(scoring,
                                                  seq_a[seq_i],
                                                  seq_b[seq_j]);

        // Long arithmetic needed to prevent overflow
        match_score[new_index]
          = (score_t)max4((long)match_score[old_index] + substitution_penalty,
                          (long)gap_a_score[old_index] + substitution_penalty,
                          (long)gap_b_score[old_index] + substitution_penalty,
                          (long)0);
      }

      if(!scoring->no_gaps)
      {
        // 2) Update gap_a_score[i][j] from position [i][j-1]
        old_index = ARR_2D_INDEX(score_width, i, j-1);

        gap_a_score[new_index]
          = (score_t)max4((long)match_score[old_index] + gap_open_penalty,
                          (long)gap_a_score[old_index] + scoring->gap_extend,
                          (long)gap_b_score[old_index] + gap_open_penalty,
                          (long)0);

        // 3) Update gap_b_score[i][j] from position [i-1][j]
        old_index = ARR_2D_INDEX(score_width, i-1, j);

        gap_b_score[new_index]
          = (score_t)max4((long)match_score[old_index] + gap_open_penalty,
                          (long)gap_a_score[old_index] + gap_open_penalty,
                          (long)gap_b_score[old_index] + scoring->gap_extend,
                          (long)0);
      }
    }
  }

  #ifdef DEBUG
  alignment_print_matrices(match_score, gap_a_score, gap_b_score,
                           length_a, length_b);
  #endif
  
  //
  // Store traceback matrix
  //
  SW_COMPUTATION* sw_computation
    = (SW_COMPUTATION*) malloc(sizeof(SW_COMPUTATION));

  sw_computation->match_score = match_score;
  sw_computation->gap_a_score = gap_a_score;
  sw_computation->gap_b_score = gap_b_score;

  sw_computation->score_width = score_width;
  sw_computation->score_height = score_height;
  
  sw_computation->scoring = scoring;
  sw_computation->seq_a = seq_a;
  sw_computation->seq_b = seq_b;

  BIT_ARRAY* bitarr = bit_array_create(arr_size);
  bit_array_set_all(bitarr);

  unsigned long *sorted_match_indices = malloc(sizeof(unsigned long)*arr_size);

  if(sorted_match_indices == NULL)
  {
    fprintf(stderr, "Error: ran out of memory.  smith_waterman_align()\n");
    fprintf(stderr, "  length(seq_a): %i; length(seq_b): %i\n",
            (int)length_a, (int)length_b);

    exit(EXIT_FAILURE);
  }

  unsigned long num_of_hits = 0;

  size_t pos;
  for(pos = 0; pos < arr_size; pos++)
  {
    if(match_score[pos] > 0)
    {
      sorted_match_indices[num_of_hits++] = pos;
    }
  }

  // Now sort matched hits
  MatrixSort tmp_struct = {match_score, score_width};
  sort_r(sorted_match_indices, num_of_hits,
         sizeof(unsigned long), sort_match_indices, &tmp_struct);

  /*
  #ifdef DEBUG
  for(i = 0; i < num_of_hits; i++)
  {
    printf(" arr_index %2lu (%li, %li): %u\n", sorted_match_indices[i],
           ARR_2D_X(sorted_match_indices[i], score_width),
           ARR_2D_Y(sorted_match_indices[i], score_width),
           match_score[sorted_match_indices[i]]);
  }
  #endif
  */

  sw_computation->match_scores_mask = bitarr;
  sw_computation->sorted_match_indices = sorted_match_indices;
  sw_computation->num_of_hits = num_of_hits;

  // Start with the first hit
  sw_computation->next_hit = 0;

  #ifdef DEBUG
  print_computation(sw_computation);
  #endif

  return sw_computation;
}

// Return 1 if alignment was found, 0 otherwise
char _follow_hit(SW_COMPUTATION* sw_computation, size_t arr_index,
                 SW_LOCAL_ALIGNMENT* result)
{
  // Follow path through matrix
  size_t score_x = (size_t)ARR_2D_X(arr_index, sw_computation->score_width);
  size_t score_y = (size_t)ARR_2D_Y(arr_index, sw_computation->score_width);

  // Local alignments always (start and) end with a match
  enum Matrix curr_matrix = MATCH;
  score_t curr_score = sw_computation->match_score[arr_index];

  // Store end arr_index and (x,y) coords for later
  size_t end_arr_index = arr_index;
  size_t end_score_x = score_x;
  size_t end_score_y = score_y;
  score_t end_score = curr_score;

  unsigned int length;

  for(length = 0; ; length++)
  {
    if(!bit_array_get_bit(sw_computation->match_scores_mask, arr_index))
    {
      //printf("  bit is zero (%i, %i)\n", score_x, score_y);
      return 0;
    }

    bit_array_clear_bit(sw_computation->match_scores_mask, arr_index);

    if(curr_score == 0)
    {
      break;
    }

    //printf(" %i (%i, %i)\n", length, score_x, score_y);

    // Find out where to go next
    alignment_reverse_move(&curr_matrix, &curr_score,
                           &score_x, &score_y, &arr_index,
                           sw_computation->score_width,
                           sw_computation->score_height,
                           sw_computation->match_score, // match matrix
                           sw_computation->gap_a_score, // gap a matrix
                           sw_computation->gap_b_score, // gap b matrix
                           sw_computation->seq_a,
                           sw_computation->seq_b,
                           sw_computation->scoring);
  }
  
  #ifdef DEBUG
  printf("assembling results...\n");
  #endif
  
  // We got a result!
  // Allocate memory for the result
  result->length = length;

  if(result->result_a == NULL)
  {
    result->result_a = (char*) malloc(length+1);
    result->result_b = (char*) malloc(length+1);
    result->capacity = length+1;
  }
  else if(result->capacity < length+1)
  {
    result->result_a = (char*) realloc(result->result_a, length+1);
    result->result_b = (char*) realloc(result->result_b, length+1);
  }

  if(result->result_a == NULL || result->result_b == NULL)
  {
    fprintf(stderr, "SmithWaterman Error: Ran out of memory for hit result\n");
    exit(-1);
  }

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
        result->result_a[i] = sw_computation->seq_a[score_x-1];
        result->result_b[i] = sw_computation->seq_b[score_y-1];
        break;

      case GAP_A:
        result->result_a[i] = '-';
        result->result_b[i] = sw_computation->seq_b[score_y-1];
        break;

      case GAP_B:
        result->result_a[i] = sw_computation->seq_a[score_x-1];
        result->result_b[i] = '-';
        break;
      
      default:
      fprintf(stderr, "Program error: invalid matrix in _follow_hit()\n");
      fprintf(stderr, "Please submit a bug report to: turner.isaac@gmail.com\n");
      exit(EXIT_FAILURE);
    }

    alignment_reverse_move(&curr_matrix, &curr_score,
                           &score_x, &score_y, &arr_index,
                           sw_computation->score_width,
                           sw_computation->score_height,
                           sw_computation->match_score, // match matrix
                           sw_computation->gap_a_score, // gap a matrix
                           sw_computation->gap_b_score, // gap b matrix
                           sw_computation->seq_a,
                           sw_computation->seq_b,
                           sw_computation->scoring);
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

// Get the next local alignment from this computation
char smith_waterman_get_hit(SW_COMPUTATION* sw_computation,
                            SW_LOCAL_ALIGNMENT* result)
{
  while(sw_computation->next_hit < sw_computation->num_of_hits)
  {
    unsigned long arr_index
      = sw_computation->sorted_match_indices[sw_computation->next_hit++];

    //printf("hit %lu/%lu\n", sw_computation->next_hit,
    //       sw_computation->num_of_hits);

    if(bit_array_get_bit(sw_computation->match_scores_mask, arr_index))
    {
      char success = _follow_hit(sw_computation, arr_index, result);

      if(success)
      {
        return 1;
      }
    }
  }

  return 0;
}

SW_LOCAL_ALIGNMENT* smith_waterman_create_hit()
{
  SW_LOCAL_ALIGNMENT* alignment
    = (SW_LOCAL_ALIGNMENT*) malloc(sizeof(SW_LOCAL_ALIGNMENT));

  alignment->result_a = NULL;
  alignment->result_b = NULL;

  return alignment;
}

void smith_waterman_free_hit(SW_LOCAL_ALIGNMENT* alignment)
{
  if(alignment->result_a != NULL)
  {
    free(alignment->result_a);
    free(alignment->result_b);
  }

  free(alignment);
}

void smith_waterman_free(SW_COMPUTATION* sw_computation)
{
  free(sw_computation->match_score);
  free(sw_computation->gap_a_score);
  free(sw_computation->gap_b_score);

  free(sw_computation->sorted_match_indices);

  bit_array_free(sw_computation->match_scores_mask);

  free(sw_computation);
}
