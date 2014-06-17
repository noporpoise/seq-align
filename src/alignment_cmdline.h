/*
 alignment_cmdline.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#ifndef ALIGNMENT_CMDLINE_HEADER_SEEN
#define ALIGNMENT_CMDLINE_HEADER_SEEN

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdarg.h> // required for va_list
#include "seq_file.h"
#include "alignment.h"

typedef struct
{
  // file inputs
  size_t file_list_length, file_list_capacity;
  char **file_paths1, **file_paths2;

  // All values initially 0
  char case_sensitive;
  int match, mismatch, gap_open, gap_extend;

  // SW specific
  score_t min_score;
  unsigned int print_context, max_hits_per_alignment;
  char min_score_set, max_hits_per_alignment_set;
  char print_seq;

  // NW specific
  char freestartgap_set, freeendgap_set;
  char print_scores;
  char zam_stle_output;

  // General output
  char print_fasta, print_pretty, print_colour;

  // Experimental
  char no_gaps_in1, no_gaps_in2;
  char no_mismatches;

  // Pair of sequences to align
  char *seq1, *seq2;
} cmdline_t;

char parse_entire_int(char *str, int *result);
char parse_entire_uint(char *str, unsigned int *result);

cmdline_t* cmdline_new(int argc, char **argv, scoring_t *scoring, char is_sw);
void cmdline_free(cmdline_t* cmd);

void cmdline_add_files(cmdline_t* cmd, char* p1, char* p2);
size_t cmdline_get_num_of_file_pairs(cmdline_t* cmd);
char* cmdline_get_file1(cmdline_t* cmd, size_t i);
char* cmdline_get_file2(cmdline_t* cmd, size_t i);

void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2));

#endif
