/*
 alignment.h
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

#ifndef ALIGNMENT_HEADER_SEEN
#define ALIGNMENT_HEADER_SEEN

#include "alignment_scoring.h"
#include "string_buffer.h"

/*
typedef struct Alignment Alignment // previously SW_LOCAL_ALIGNMENT
*/

// C Preprocessor sets SCORE_TYPE
// unsigned int for SW; int for NW
typedef SCORE_TYPE score_t;

#define MATRIX_NAME(x) ((x) == MATCH ? "MATCH" : ((x) == GAP_A ? "GAP_A" : "GAP_B"))

// Matrix names
enum Matrix { MATCH,GAP_A,GAP_B };

// Printing colour codes
char *align_col_mismatch, *align_col_indel, *align_col_context, *align_col_stop;

/*
// Constructors/Destructors
Alignment* alignment_create(size_t capacity);
void alignment_ensure_capacity(Alignment* alignment);
void alignment_destroy(Alignment* alignment);
*/

long max2(long a, long b);
long max3(long a, long b, long c);
long max4(long a, long b, long c, long d);

// Methods
void alignment_print_matrices(const score_t* match_score,
                              const score_t* gap_a_score,
                              const score_t* gap_b_score,
                              int length_a, int length_b);

void alignment_colour_print_against(const char *alignment_a,
                                    const char *alignment_b,
                                    char case_sensitive);

void alignment_print_spacer(const char* alignment_a, const char* alignment_b,
                            const SCORING_SYSTEM* scoring);

void alignment_reverse_move(enum Matrix *curr_matrix, score_t* curr_score,
                            size_t *score_x, size_t *score_y,
                            unsigned long *arr_index,
                            size_t score_width, size_t score_height,
                            const score_t *match_score,
                            const score_t *gap_a_score,
                            const score_t *gap_b_score,
                            const char* seq_a, const char* seq_b,
                            const SCORING_SYSTEM* scoring);

#endif
