/*
 alignment_scoring_load.h
 project: AlignmentScoring
 author: Isaac Turner <turner.isaac@gmail.com>
 Used in SmithWaterman and NeedlemanWunsch projects
 url: http://sourceforge.net/projects/needlemanwunsch
 url: http://sourceforge.net/projects/smithwaterman
 Copyright (C) 06-Dec-2011
 
 see: README

 == License
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

#ifndef ALIGNMENT_SCORING_LOAD_HEADER_SEEN
#define ALIGNMENT_SCORING_LOAD_HEADER_SEEN

#include "alignment_scoring.h"

void align_scoring_load_matrix(gzFile* file, const char* file_path,
                               SCORING_SYSTEM* scoring, char case_sensitive);

void align_scoring_load_pairwise(gzFile* file, const char* file_path,
                                 SCORING_SYSTEM* scoring, char case_sensitive);

#endif
