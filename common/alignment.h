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

#define MATRIX_NAME(x) ((x) == MATCH ? "MATCH" : ((x) == GAP_A ? "GAP_A" : "GAP_B"))

typedef SCORE_TYPE score_t;

// Matrix names
enum Matrix { MATCH,GAP_A,GAP_B };

// Printing colour codes
char *align_col_mismatch, *align_col_indel, *align_col_context, *align_col_stop;

// Methods
void alignment_print_matrices(score_t* match_score, score_t* gap_a_score,
                              score_t* gap_b_score,
                              int length_a, int length_b);

void alignment_colour_print_against(const char *alignment_a,
                                    const char *alignment_b,
                                    char case_sensitive);

void alignment_print_spacer(const char* alignment_a, const char* alignment_b,
                            const SCORING_SYSTEM* scoring);

void alignment_reverse_move(enum Matrix *curr_matrix, score_t* curr_score,
                            unsigned int *score_x, unsigned int *score_y,
                            unsigned long *arr_index,
                            unsigned int score_width,
                            const score_t *match_score,
                            const score_t *gap_a_score,
                            const score_t *gap_b_score,
                            const char* seq_a, const char* seq_b,
                            const SCORING_SYSTEM* scoring);

#endif
