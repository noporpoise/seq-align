/*
 tools/lcs_cmdline.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Feb 2015
 */

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower

#include "smith_waterman.h"

static void print_usage(char **argv)
{
  fprintf(stderr, "%s [options] <sequence>\n", argv[0]);
  fprintf(stderr, "  Print substrings in decreasing order of length\n");
//   fprintf(stderr, "Options:\n"
// "  -c,--case-insensitive  Case insensitive matching\n"
// "  -C,--case-sensitive    Case sensitive matching [default]\n"
// "  -p,--positions         Print positions\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
  if(argc != 2) print_usage(argv);

  char *seq = argv[1];
  size_t seqlen = strlen(seq);

  // Go
  int match = 1, mismatch = -1, gap_open = -4, gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = true, no_gaps_in_b = true;
  bool no_mismatches = true, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  // Alignment results stored here
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *aln = alignment_create(seqlen+1);

  smith_waterman_align(seq, seq, &scoring, sw);

  // Loop over results
  while(smith_waterman_fetch(sw, aln))
  {
    if(aln->pos_a < aln->pos_b) {
      fputs(aln->result_a, stdout);
      printf(" [%zu,%zu]\n", aln->pos_a, aln->pos_b);
    }
  }

  smith_waterman_free(sw);
  alignment_free(aln);

  return EXIT_SUCCESS;
}
