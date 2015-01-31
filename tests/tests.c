#include "cutest.h"
#include "needleman_wunsch.h"
#include <string.h>

//No gap is expected in the longer sequence
void test_no_gaps_in_longer(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "aaaaacg";
  const char* seq_b = "acgt"; 

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  
  char no_gaps_in_a = 1, no_gaps_in_b = 0;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               0, 0,
               no_gaps_in_a, no_gaps_in_b, 0, 0, 0, 0);
               
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_a, "aaaaacg") == 0 && strcmp(result->result_b, "acgt---")== 0);
  
 }
 
 //If the sequences have same length "n", with nogaps the length of alignment is "n"
 void test_no_gaps_equal(void)
 {
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "actcg";
  const char* seq_b = "ggctc"; 

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  
  char no_gaps_in_a = 1, no_gaps_in_b = 1;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               0, 0,
               no_gaps_in_a, no_gaps_in_b, 0, 0, 0, 0);           
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_a, "actcg") == 0 && strcmp(result->result_b, "ggctc")== 0);
}

void test_gaps_only_at_ends_in1(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "acgtc";
  const char* seq_b = "aaaaacgc"; 

  // Decide on scoring
  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  
  char gaps_ends_a = 1, no_gaps_in_a = 0, no_gaps_in_b = 0;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               0, 0,
               no_gaps_in_a, no_gaps_in_b, gaps_ends_a, 0, 0, 0);
               
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_b, "aaaaacg-c") == 0 && strcmp(result->result_a, "----acgtc")== 0);
  seq_b = "aaaaacgct";
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_b, "aaaaacgct") == 0 && strcmp(result->result_a, "----acgtc")== 0);
  mismatch = -10;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               0, 0,
               no_gaps_in_a, no_gaps_in_b, gaps_ends_a, 0, 0, 0);
               
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_b, "aaaaacg-ct") == 0 && strcmp(result->result_a, "----acgtc-")== 0);
}

TEST_LIST = {
  { "no_gaps_in_longer", test_no_gaps_in_longer },
  { "no_gaps_equal_length", test_no_gaps_equal},
  { "gaps_only_at_ends_in_shorter", test_gaps_only_at_ends_in1},
  { 0 }
};