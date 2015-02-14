/*
 examples/nw_example.cpp
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Feb 2015
 */


#include <iostream>
#include "needleman_wunsch.h"

int main(int argc, char* argv[])
{
  using namespace std;

  if(argc != 3)
  {
    cout << "Usage: " << argv[0] << " <seq1> <seq2>\n";
    return -1;
  }

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
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);

  needleman_wunsch_align(argv[1], argv[2], &scoring, nw, aln);

  cout << "seqA: " << aln->result_a << "\n";
  cout << "seqB: " << aln->result_b << "\n";
  cout << "alignment score: " << aln->score << "\n";

  needleman_wunsch_free(nw);
  alignment_free(aln);

  return 0;
}
