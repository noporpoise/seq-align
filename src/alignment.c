/*
 alignment.c
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 Nov 2013; License: Public Domain
 */

// Turn on debugging output by defining DEBUG
//#define DEBUG

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h> // INT_MIN
#include <ctype.h> // tolower

#include "alignment.h"
#include "alignment_macros.h"

const char align_col_mismatch[] = "\033[92m"; // Mismatch (GREEN)
const char align_col_indel[] = "\033[91m"; // Insertion / deletion (RED)
// Pink used by SmithWaterman local alignment for printing surrounding bases
const char align_col_context[] = "\033[95m";
const char align_col_stop[] = "\033[0m";

static long max4(long a, long b, long c, long d)
{
  long result = a;
  if(b > result) result = b;
  if(c > result) result = c;
  if(d > result) result = d;
  return result;
}

// Fill in traceback matrix
static void alignment_fill_matrices(aligner_t *aligner, char is_sw)
{
  score_t *match_scores = aligner->match_scores;
  score_t *gap_a_scores = aligner->gap_a_scores;
  score_t *gap_b_scores = aligner->gap_b_scores;
  const scoring_t *scoring = aligner->scoring;
  size_t score_width = aligner->score_width;
  size_t score_height = aligner->score_height;

  size_t i, j, arr_size = score_width * score_height;

  score_t min = is_sw ? 0 : INT_MIN;

  size_t seq_i, seq_j, len_i = score_width-1, len_j = score_height-1;
  size_t index, index_left, index_up, index_upleft;

  if(scoring->no_gaps_in_a || scoring->gaps_ends_a) {
    for(i = 0; i < arr_size; i++) gap_a_scores[i] = min;
  }

  if(scoring->no_gaps_in_b || scoring->gaps_ends_b) {
    for(i = 0; i < arr_size; i++) gap_b_scores[i] = min;
  }

  // [0][0]
  match_scores[0] = 0;
  gap_a_scores[0] = 0;
  gap_b_scores[0] = 0;

  // work along first row -> [i][0]
  for(i = 1; i < score_width; i++)
  {
    match_scores[i] = min;

    // Think carefully about which way round these are
    gap_a_scores[i] = min;
    if((len_j < len_i && scoring->no_gaps_in_b) || !scoring->no_gaps_in_b)
   	 gap_b_scores[i] = (score_t)(scoring->no_start_gap_penalty ? 0
                      : scoring->gap_open + (long)i * scoring->gap_extend);
  }

  // work down first column -> [0][j]
  for(j = 1, index = score_width; j < score_height; j++, index += score_width)
  {
    match_scores[index] = min;

    // Think carefully about which way round these are
    if((len_i < len_j && scoring->no_gaps_in_a) || !scoring->no_gaps_in_a)
    	gap_a_scores[index]
      		= (score_t)(scoring->no_start_gap_penalty ? 0
                  : scoring->gap_open + (long)j * scoring->gap_extend);
    gap_b_scores[index] = min;
  }

  // These are longs to force addition to be done with higher accuracy
  long gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  long gap_extend_penalty = scoring->gap_extend;
  long substitution_penalty;

  // start at position [1][1]
  index_upleft = 0;
  index_up = 1;
  index_left = score_width;
  index = score_width+1;

  for(seq_j = 0; seq_j < len_j; seq_j++)
  {
    for(seq_i = 0; seq_i < len_i; seq_i++)
    {
      // Update match_scores[i][j] with position [i-1][j-1]
      // substitution penalty
      char is_match;
      int tmp_penalty;

      scoring_lookup(scoring, aligner->seq_a[seq_i], aligner->seq_b[seq_j],
                     &tmp_penalty, &is_match);

      if(scoring->no_mismatches && !is_match)
      {
        match_scores[index] = min;
      }
      else
      {
        substitution_penalty = tmp_penalty; // cast to long

        // substitution
        // 1) continue alignment
        // 2) close gap in seq_a
        // 3) close gap in seq_b
        match_scores[index]
          = max4(match_scores[index_upleft] + substitution_penalty,
                 gap_a_scores[index_upleft] + substitution_penalty,
                 gap_b_scores[index_upleft] + substitution_penalty,
                 min);
      }

      // Long arithmetic since some INTs are set to min and penalty is -ve
      // (adding as ints would cause an integer overflow)

      if(!scoring->no_gaps_in_a && !scoring->gaps_ends_a)
      {
        // Update gap_a_scores[i][j] from position [i][j-1]

        if(seq_i == len_i-1 && scoring->no_end_gap_penalty)
        {
          gap_a_scores[index] = MAX3(match_scores[index_up],
                                     gap_a_scores[index_up],
                                     gap_b_scores[index_up]);
        }
        else
        {
          gap_a_scores[index]
            = max4(match_scores[index_up] + gap_open_penalty,
                   gap_a_scores[index_up] + gap_extend_penalty,
                   gap_b_scores[index_up] + gap_open_penalty,
                   min);
        }
      }

      if(!scoring->no_gaps_in_b && !scoring->gaps_ends_b)
      {
        // Update gap_b_scores[i][j] from position [i-1][j]

        if(seq_j == len_j-1 && scoring->no_end_gap_penalty)
        {
          gap_b_scores[index] = MAX3(match_scores[index_left],
                                     gap_a_scores[index_left],
                                     gap_b_scores[index_left]);
        }
        else
        {
          gap_b_scores[index]
            = max4(match_scores[index_left] + gap_open_penalty,
                   gap_a_scores[index_left] + gap_open_penalty,
                   gap_b_scores[index_left] + gap_extend_penalty,
                   min);
        }
      }

      index++;
      index_left++;
      index_up++;
      index_upleft++;
    }

    index++;
    index_left++;
    index_up++;
    index_upleft++;
  }

  if((scoring->no_gaps_in_a && len_i < len_j) || scoring->gaps_ends_a)
  {
    // Allow gaps only at the start/end of A
    // Fill right hand column of traceback matrix
    index_up = score_width-1;
    index = index_up+score_width;

    while(index < arr_size)
    {
      gap_a_scores[index]
        = MAX3(match_scores[index_up] + ((scoring->no_end_gap_penalty) ? 0 : gap_open_penalty),
               gap_a_scores[index_up] + ((scoring->no_end_gap_penalty) ? 0 : gap_extend_penalty),
               min);

      index_up = index;
      index += score_width;
    }
  }

  if((scoring->no_gaps_in_b && len_j < len_i) || scoring->gaps_ends_b)
  {
    // Allow gaps only at the start/end of B
    // Fill bottom row of traceback matrix
    index_left = (score_height-1)*score_width;
    index = index_left+1;

    while(index < arr_size)
    {
      gap_b_scores[index]
        = MAX3(match_scores[index_left] + ((scoring->no_end_gap_penalty) ? 0 : gap_open_penalty),
               gap_b_scores[index_left] + ((scoring->no_end_gap_penalty) ? 0 : gap_extend_penalty),
               min);

      index_left++;
      index++;
    }
  }
}

void aligner_align(aligner_t *aligner,
                   const char *seq_a, const char *seq_b,
                   size_t len_a, size_t len_b,
                   const scoring_t *scoring, char is_sw)
{
  aligner->scoring = scoring;
  aligner->seq_a = seq_a;
  aligner->seq_b = seq_b;
  aligner->score_width = len_a+1;
  aligner->score_height = len_b+1;

  size_t new_capacity = aligner->score_width * aligner->score_height;

  if(aligner->capacity == 0 || aligner->capacity < new_capacity)
  {
    aligner->capacity = ROUNDUP2POW(new_capacity);
    size_t mem = sizeof(score_t) * aligner->capacity;

    if(aligner->capacity == 0) {
      aligner->match_scores = malloc(mem);
      aligner->gap_a_scores = malloc(mem);
      aligner->gap_b_scores = malloc(mem);
    } else {
      aligner->match_scores = realloc(aligner->match_scores, mem);
      aligner->gap_a_scores = realloc(aligner->gap_a_scores, mem);
      aligner->gap_b_scores = realloc(aligner->gap_b_scores, mem);
    }
  }

  alignment_fill_matrices(aligner, is_sw);
}

void aligner_destroy(aligner_t *aligner)
{
  if(aligner->capacity > 0) {
    free(aligner->match_scores);
    free(aligner->gap_a_scores);
    free(aligner->gap_b_scores);
  }
}


alignment_t* alignment_create(size_t capacity)
{
  capacity = ROUNDUP2POW(capacity);
  alignment_t *result = malloc(sizeof(alignment_t));
  result->result_a = malloc(sizeof(char)*capacity);
  result->result_b = malloc(sizeof(char)*capacity);
  result->capacity = capacity;
  result->length = 0;
  result->result_a[0] = result->result_b[0] = '\0';
  result->pos_a = result->pos_b = result->len_a = result->len_b = 0;
  result->score = 0;
  return result;
}

void alignment_ensure_capacity(alignment_t* result, size_t strlength)
{
  size_t capacity = strlength+1;
  if(result->capacity < capacity)
  {
    capacity = ROUNDUP2POW(capacity);
    result->result_a = realloc(result->result_a, sizeof(char)*capacity);
    result->result_b = realloc(result->result_b, sizeof(char)*capacity);
    result->capacity = capacity;
    if(result->result_a == NULL || result->result_b == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

void alignment_free(alignment_t* result)
{
  free(result->result_a);
  free(result->result_b);
}


// Backtrack through scoring matrices
void alignment_reverse_move(enum Matrix *curr_matrix, score_t *curr_score,
                            size_t *score_x, size_t *score_y,
                            size_t *arr_index, const aligner_t *aligner)
{
  size_t seq_x = (*score_x)-1, seq_y = (*score_y)-1;
  size_t len_i = aligner->score_width-1, len_j = aligner->score_height-1;

  char is_match;
  int match_penalty;
  const scoring_t *scoring = aligner->scoring;

  scoring_lookup(scoring, aligner->seq_a[seq_x], aligner->seq_b[seq_y],
                 &match_penalty, &is_match);

  int gap_a_open_penalty, gap_b_open_penalty;
  int gap_a_extend_penalty, gap_b_extend_penalty;

  gap_a_open_penalty = gap_b_open_penalty = scoring->gap_extend + scoring->gap_open;
  gap_a_extend_penalty = gap_b_extend_penalty = scoring->gap_extend;

  // Free gaps at the ends
  if(scoring->no_end_gap_penalty) {
    if(*score_x == len_i) gap_a_open_penalty = gap_a_extend_penalty = 0;
    if(*score_y == len_j) gap_b_open_penalty = gap_b_extend_penalty = 0;
  }
  if(scoring->no_start_gap_penalty) {
    if(*score_x == 0) gap_a_open_penalty = gap_a_extend_penalty = 0;
    if(*score_y == 0) gap_b_open_penalty = gap_b_extend_penalty = 0;
  }

  long prev_match_penalty, prev_gap_a_penalty, prev_gap_b_penalty;

  switch(*curr_matrix)
  {
    case MATCH:
      prev_match_penalty = match_penalty;
      prev_gap_a_penalty = match_penalty;
      prev_gap_b_penalty = match_penalty;
      (*score_x)--;
      (*score_y)--;
      (*arr_index) -= aligner->score_width + 1;
      break;

    case GAP_A:
      prev_match_penalty = gap_a_open_penalty;
      prev_gap_a_penalty = gap_a_extend_penalty;
      prev_gap_b_penalty = gap_a_open_penalty;
      (*score_y)--;
      (*arr_index) -= aligner->score_width;
      break;

    case GAP_B:
      prev_match_penalty = gap_b_open_penalty;
      prev_gap_a_penalty = gap_b_open_penalty;
      prev_gap_b_penalty = gap_b_extend_penalty;
      (*score_x)--;
      (*arr_index)--;
      break;

    default:
      fprintf(stderr, "Program error: invalid matrix in get_reverse_move()\n");
      fprintf(stderr, "Please submit a bug report to: turner.isaac@gmail.com\n");
      exit(EXIT_FAILURE);
  }

  // *arr_index = ARR_2D_INDEX(aligner->score_width, *score_x, *score_y);

  if((!scoring->no_gaps_in_a || *score_x == 0 || *score_x == len_i) &&
     (long)aligner->gap_a_scores[*arr_index] + prev_gap_a_penalty == *curr_score)
  {
    *curr_matrix = GAP_A;
    *curr_score = aligner->gap_a_scores[*arr_index];
  }
  else if((!scoring->no_gaps_in_b || *score_y == 0 || *score_y == len_j) &&
          (long)aligner->gap_b_scores[*arr_index] + prev_gap_b_penalty == *curr_score)
  {
    *curr_matrix = GAP_B;
    *curr_score = aligner->gap_b_scores[*arr_index];
  }
  else if((!scoring->no_mismatches || is_match) || scoring->no_mismatches &&
          (long)aligner->match_scores[*arr_index] + prev_match_penalty == *curr_score)
  {
    *curr_matrix = MATCH;
    *curr_score = aligner->match_scores[*arr_index];
  }
  else
  {
    fprintf(stderr, "[%s:%zu,%zu]: %i\n",
            MATRIX_NAME(*curr_matrix), *score_x, *score_y, *curr_score);
    fprintf(stderr, " match: %li gap_open: %li gap_extend: %li\n",
            prev_match_penalty, prev_gap_a_penalty, prev_gap_b_penalty);
    fprintf(stderr, " MATCH: %i GAP_A: %i GAP_B: %i\n",
            aligner->match_scores[*arr_index],
            aligner->gap_a_scores[*arr_index],
            aligner->gap_b_scores[*arr_index]);

    fprintf(stderr,
"Program error: traceback fail (get_reverse_move)\n"
"This may be due to an integer overflow if your sequences are long or scores\n"
"are large. If this is the case using smaller scores or shorter sequences may\n"
"work around this problem.  \n"
"  If you think this is a bug, please report it to: turner.isaac@gmail.com\n");
    exit(EXIT_FAILURE);
  }
}


void alignment_print_matrices(const aligner_t *aligner)
{
  const score_t* match_scores = aligner->match_scores;
  const score_t* gap_a_scores = aligner->gap_a_scores;
  const score_t* gap_b_scores = aligner->gap_b_scores;

  size_t i, j;

  printf("match_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf(" %3i", (int)ARR_LOOKUP(match_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }
  printf("gap_a_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf(" %3i", (int)ARR_LOOKUP(gap_a_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }
  printf("gap_b_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf(" %3i", (int)ARR_LOOKUP(gap_b_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
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
        fputs(align_col_indel, stdout);
        red = 1;
      }
    }
    else if(red)
    {
      red = 0;
      fputs(align_col_stop, stdout);
    }

    if(((case_sensitive && alignment_a[i] != alignment_b[i]) ||
        (!case_sensitive && tolower(alignment_a[i]) != tolower(alignment_b[i]))) &&
       alignment_a[i] != '-' && alignment_b[i] != '-')
    {
      if(!green)
      {
        fputs(align_col_mismatch, stdout);
        green = 1;
      }
    }
    else if(green)
    {
      green = 0;
      fputs(align_col_stop, stdout);
    }

    putc(alignment_a[i], stdout);
  }

  if(green || red)
  {
    // Stop all colours
    fputs(align_col_stop, stdout);
  }
}

// Order of alignment_a / alignment_b is not important
void alignment_print_spacer(const char* alignment_a, const char* alignment_b,
                            const scoring_t* scoring)
{
  int i;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-' || alignment_b[i] == '-')
    {
      putc(' ', stdout);
    }
    else if(alignment_a[i] == alignment_b[i] ||
            (!scoring->case_sensitive &&
             tolower(alignment_a[i]) == tolower(alignment_b[i])))
    {
      putc('|', stdout);
    }
    else
    {
      putc('*', stdout);
    }
  }
}
