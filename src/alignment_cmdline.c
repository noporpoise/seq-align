/*
 alignment_cmdline.c
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
#include <limits.h>

#include "alignment_cmdline.h"

#include "seq_file.h"

// File loading
int file_list_length = 0;
int file_list_capacity = 0;
char **file_paths1 = NULL, **file_paths2 = NULL;

char parse_entire_int(char *str, int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (int)tmp;
    return 1;
  }
}

char parse_entire_uint(char *str, unsigned int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(tmp > UINT_MAX || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (unsigned int)tmp;
    return 1;
  }
}

void cmdline_init()
{
  file_paths1 = (char**)malloc(sizeof(char*)*file_list_capacity);
  file_paths2 = (char**)malloc(sizeof(char*)*file_list_capacity);
}

void cmdline_finish()
{
  free(file_paths1);
  free(file_paths2);
}

void _check_file_array_lengths()
{
  if(file_list_capacity == 0)
  {
    file_list_capacity = 10;
    file_paths1 = malloc(sizeof(char*) * file_list_capacity);
    file_paths2 = malloc(sizeof(char*) * file_list_capacity);

    if(file_paths1 == NULL || file_paths2 == NULL)
    {
      fprintf(stderr, "%s:%i: Ran out of memory taking file arguments!\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
  else if(file_list_length == file_list_capacity)
  {
    // Expand arrays used for holding file paths
    file_list_capacity *= 2;
    file_paths1 = realloc(file_paths1, sizeof(char*)*file_list_capacity);
    file_paths2 = realloc(file_paths2, sizeof(char*)*file_list_capacity);

    if(file_paths1 == NULL || file_paths2 == NULL)
    {
      fprintf(stderr, "%s:%i: Ran out of memory taking file arguments!\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

void cmdline_add_files(char* p1, char* p2)
{
  _check_file_array_lengths();

  file_paths1[file_list_length] = p1;
  file_paths2[file_list_length] = p2;

  file_list_length++;
}

int cmdline_get_num_of_file_pairs()
{
  return file_list_length;
}

char* cmdline_get_file1(int i)
{
  return file_paths1[i];
}

char* cmdline_get_file2(int i)
{
  return file_paths2[i];
}

// If seq2 is NULL, read pair of entries from first file
// Otherwise read an entry from each
void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2))
{
  seq_file_t *sf1, *sf2;

  if((sf1 = seq_open(path1)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path1);
    fflush(stderr);
    return;
  }

  if(path2 == NULL)
  {
    sf2 = sf1;
  }
  else if((sf2 = seq_open(path2)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path1);
    fflush(stderr);
    return;
  }

  read_t *read1 = seq_read_alloc();
  read_t *read2 = seq_read_alloc();

  // Loop while we can read a sequence from the first file
  unsigned long alignments;

  for(alignments = 0; seq_read(sf1, read1) > 0; alignments++)
  {
    if(seq_read(sf2, read2) <= 0)
    {
      fprintf(stderr, "Alignment Error: Odd number of sequences - "
                      "I read in pairs!\n");
      fflush(stderr);
      break;
    }

    (align)(read1, read2);
  }

  // warn if no bases read
  if(alignments == 0)
  {
    fprintf(stderr, "Alignment Warning: empty input\n");
    fflush(stderr);
  }

  // Close files
  seq_close(sf1);

  if(path2 != NULL)
    seq_close(sf2);

  // Free memory
  seq_read_destroy(read1);
  seq_read_destroy(read2);
}
