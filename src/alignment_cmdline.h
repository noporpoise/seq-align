/*
 alignment_cmdline.h
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

#ifndef ALIGNMENT_CMDLINE_HEADER_SEEN
#define ALIGNMENT_CMDLINE_HEADER_SEEN

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdarg.h> // required for va_list
#include <stdbool.h>
#include "seq_file/seq_file.h"
#include "alignment.h"

enum SeqAlignCmdType {SEQ_ALIGN_SW_CMD, SEQ_ALIGN_NW_CMD, SEQ_ALIGN_LCS_CMD};

typedef struct
{
  // file inputs
  size_t file_list_length, file_list_capacity;
  char **file_paths1, **file_paths2;

  // All values initially 0
  bool case_sensitive;
  int match, mismatch, gap_open, gap_extend;

  // SW specific
  score_t min_score;
  unsigned int print_context, max_hits_per_alignment;
  bool min_score_set, max_hits_per_alignment_set;
  bool print_seq;

  // NW specific
  bool freestartgap_set, freeendgap_set;
  bool print_matrices, print_scores;
  bool zam_stle_output;

  // Turns off zlib for stdin
  bool interactive;

  // General output
  bool print_fasta, print_pretty, print_colour;

  // Experimental
  bool no_gaps_in1, no_gaps_in2;
  bool no_mismatches;

  // Pair of sequences to align
  const char *seq1, *seq2;
} cmdline_t;

char parse_entire_int(char *str, int *result);
char parse_entire_uint(char *str, unsigned int *result);

cmdline_t* cmdline_new(int argc, char **argv, scoring_t *scoring,
                       enum SeqAlignCmdType cmd_type);
void cmdline_free(cmdline_t* cmd);

void cmdline_add_files(cmdline_t* cmd, char* p1, char* p2);
size_t cmdline_get_num_of_file_pairs(cmdline_t* cmd);
char* cmdline_get_file1(cmdline_t* cmd, size_t i);
char* cmdline_get_file2(cmdline_t* cmd, size_t i);

void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2),
                     bool use_zlib);

#endif
