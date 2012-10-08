/*
 smith_waterman.h
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

#ifndef SMITH_WATERMAN_HEADER_SEEN
#define SMITH_WATERMAN_HEADER_SEEN

#include "bit_array.h"
#include "alignment.h"

/*
typedef struct SWAligner      SWAligner
typedef struct SWAlignment    SWAlignment    // previously SW_COMPUTATION
*/

typedef struct SW_COMPUTATION SW_COMPUTATION;
typedef struct SW_LOCAL_ALIGNMENT SW_LOCAL_ALIGNMENT;

struct SW_LOCAL_ALIGNMENT
{
  // Store local alignment result here
  char *result_a, *result_b;
  unsigned int capacity, length;
  size_t pos_a, pos_b; // position of first base (0-based)
  size_t len_a, len_b; // number of bases in alignment
  score_t score;
};

/*
 Do not alter seq_a, seq_b or scoring whilst calling this method
 or between calls to smith_waterman_get_hit
*/
SW_COMPUTATION* smith_waterman_align(const char* seq_a, const char* seq_b,
                                     SCORING_SYSTEM* scoring);

size_t smith_waterman_seq_a_strlen(SW_COMPUTATION *sw);
size_t smith_waterman_seq_b_strlen(SW_COMPUTATION *sw);

SW_LOCAL_ALIGNMENT* smith_waterman_create_hit();

// An alignment to read from, and a pointer to memory to store the result
// returns 1 if an alignment was read, 0 otherwise
char smith_waterman_get_hit(SW_COMPUTATION* sw_computation,
                            SW_LOCAL_ALIGNMENT* result);

void smith_waterman_free_hit(SW_LOCAL_ALIGNMENT* alignment);
void smith_waterman_free(SW_COMPUTATION* sw_computation);

#endif
