/*
 needleman_wunsch.h
 project: NeedlemanWunsch
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/needlemanwunsch
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


#ifndef NEEDLEMAN_WUNSCH_HEADER_SEEN
#define NEEDLEMAN_WUNSCH_HEADER_SEEN

#include "alignment.h"

// alloc memory for result (returns length of seq_a + seq_b)
int nw_alloc_mem(const char* seq_a, const char* seq_b,
                 char** alignment_a, char** alignment_b);

// length is = length_a + length_b
// Returns 1 on success, 0 on failure
char nw_realloc_mem(unsigned int length, char** alignment_a, char** alignment_b);

/* Alignment */

int needleman_wunsch(const char* seq_a, const char* seq_b,
                     char* alignment_a, char* alignment_b,
                     SCORING_SYSTEM* scoring);

#endif
