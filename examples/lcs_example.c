#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "needleman_wunsch.h"
#include "print_lcs.h"

void align(char* seq_a, char* seq_b)
{
  // Variables to store alignment result
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  scoring_t scoring;
  scoring_system_lcs(&scoring);

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);

  print_lcs(result);

  // Free memory for storing alignment results
  needleman_wunsch_free(nw);
  alignment_free(result);
}

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    printf("usage: ./nw_example <seq1> <seq2>\n");
    exit(EXIT_FAILURE);
  }

  align(argv[1], argv[2]);
  exit(EXIT_SUCCESS);
}
