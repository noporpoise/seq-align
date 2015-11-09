/*
 tools/nw_cmdline.h
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower

// Alignment scoring and loading
#include "alignment_cmdline.h"

#include "needleman_wunsch.h"

// For this run
cmdline_t *cmd;
scoring_t scoring;

// Alignment results stored here
nw_aligner_t *nw;
alignment_t *result;

static void nw_set_default_scoring()
{
  scoring_system_default(&scoring);
}

static void align_zam(const char *seq_a, const char *seq_b)
{
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);

  // Swap '-' for '_'
  int i;
  for(i = 0; result->result_a[i] != '\0'; i++)
  {
    if(result->result_a[i] == '-') result->result_a[i] = '_';
    if(result->result_b[i] == '-') result->result_b[i] = '_';
  }

  int num_of_mismatches = 0;
  int num_of_indels = 0;

  // Print branch 1 and spacer
  printf("Br1:%s\n    ", result->result_a);

  for(i = 0; result->result_a[i] != '\0'; i++)
  {
    if(result->result_a[i] == '_' || result->result_b[i] == '_')
    {
      putc(' ', stdout);
      num_of_indels++;
    }
    else if((scoring.case_sensitive && result->result_a[i] != result->result_b[i]) ||
            tolower(result->result_a[i]) != tolower(result->result_b[i]))
    {
      putc('*', stdout);
      num_of_mismatches++;
    }
    else
    {
      putc('|', stdout);
    }
  }

  // Print branch 2 and mismatch indel numbers
  printf("\nBr2:%s\n%i %i\n\n", result->result_b,
         num_of_mismatches, num_of_indels);
}

static void align(const char *seq_a, const char *seq_b,
                  const char *seq_a_name, const char *seq_b_name)
{
  if(cmd->zam_stle_output)
  {
    align_zam(seq_a, seq_b);
    fflush(stdout);
    return;
  }

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);

  if(cmd->print_matrices)
  {
    alignment_print_matrices(nw);
  }

  if(cmd->print_fasta && seq_a_name != NULL)
  {
    fputs(seq_a_name, stdout);
    putc('\n', stdout);
  }

  if(cmd->print_fasta && cmd->print_pretty && seq_b_name != NULL)
  {
    fputs(seq_b_name, stdout);
    putc('\n', stdout);
  }

  if(cmd->print_colour)
  {
    // Print alignment line 1
    alignment_colour_print_against(result->result_a, result->result_b,
                                   scoring.case_sensitive);
  }
  else
  {
    fputs(result->result_a, stdout);
  }
  putc('\n', stdout);

  if(cmd->print_pretty)
  {
    // Print spacer
    alignment_print_spacer(result->result_a, result->result_b, &scoring);
    putc('\n', stdout);
  }
  else if(cmd->print_fasta && seq_b_name != NULL)
  {
    fputs(seq_b_name, stdout);
    putc('\n', stdout);
  }

  if(cmd->print_colour)
  {
    // Print alignment line 2
    alignment_colour_print_against(result->result_b, result->result_a,
                                   scoring.case_sensitive);
  }
  else
  {
    fputs(result->result_b, stdout);
  }
  putc('\n', stdout);

  if(cmd->print_scores) {
    printf("score: %i\n", result->score);
  }

  putc('\n', stdout);
  fflush(stdout);
}

static void align_pair_from_file(read_t *read1, read_t *read2)
{
  align(read1->seq.b, read2->seq.b,
        (read1->name.end == 0 ? NULL : read1->name.b),
        (read2->name.end == 0 ? NULL : read2->name.b));
}

int main(int argc, char* argv[])
{
  #ifdef SEQ_ALIGN_VERBOSE
  printf("VERBOSE: on\n");
  #endif

  nw_set_default_scoring();
  cmd = cmdline_new(argc, argv, &scoring, SEQ_ALIGN_NW_CMD);

  // Align!
  nw = needleman_wunsch_new();
  result = alignment_create(256);

  if(cmd->seq1 != NULL)
  {
    // Align seq1 and seq2 pair passed on the command line
    align(cmd->seq1, cmd->seq2, NULL, NULL);
  }

  // Align from files
  size_t i, num_of_file_pairs = cmdline_get_num_of_file_pairs(cmd);
  for(i = 0; i < num_of_file_pairs; i++)
  {
    const char *file1 = cmdline_get_file1(cmd, i);
    const char *file2 = cmdline_get_file2(cmd, i);
    if(file1 != NULL && *file1 == '\0' && file2 == NULL) {
      file1 = "-";
    }
    align_from_file(file1, file2, &align_pair_from_file, !cmd->interactive);
  }

  // Free memory for storing alignment results
  needleman_wunsch_free(nw);
  alignment_free(result);

  cmdline_free(cmd);

  return EXIT_SUCCESS;
}
