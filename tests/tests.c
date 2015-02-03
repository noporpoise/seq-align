#include "cutest.h"
#include "needleman_wunsch.h"
#include "smith_waterman.h"
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

/* If sequence with only gaps at ends allowed is the shorter, then the expected result is the
 * same as nogaps */
void test_gaps_only_at_ends_in_shorter(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "acgtc";
  const char* seq_b = "aaaaacgc"; 

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
  TEST_CHECK(strcmp(result->result_a, "----acgtc") == 0 && strcmp(result->result_b, "aaaaacg-c")== 0);
}

/* If sequence with only gaps at ends allowed is the longer, then the expected result is not
 * the same as nogaps because gaps are allowed at ends of the longest sequence*/
void test_gaps_only_at_ends_in_longer(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "acgct";
  const char* seq_b = "aaaaacgc"; 

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  
  char gaps_ends_b = 1, no_gaps_in_a = 0, no_gaps_in_b = 0;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               0, 0,
               no_gaps_in_a, no_gaps_in_b, 0, gaps_ends_b, 0, 0);
               
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_a, "a----cgct") == 0 && strcmp(result->result_b, "aaaaacgc-")== 0);
}

void test_no_gaps_smith_waterman(void)
{
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *result = alignment_create(256);

  const char* seq_a = "gacag";
  const char* seq_b = "tgaagt"; 
  
  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  
  char no_gaps_in_a = 1, no_gaps_in_b = 1;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               0, 0,
               no_gaps_in_a, no_gaps_in_b, 0, 0, 0, 0);
  smith_waterman_align(seq_a, seq_b, &scoring, sw);
  smith_waterman_fetch(sw, result);
  TEST_CHECK(strcmp(result->result_a, "ga") == 0 && strcmp(result->result_b, "ga") == 0);
  smith_waterman_fetch(sw, result);
  TEST_CHECK(strcmp(result->result_a, "ag") == 0 && strcmp(result->result_b, "ag") == 0);
}

/* First sequence is aligned to the end of the second one because of the reduced cost of
 * mismatches and cost-free start gap */
void test_free_start_gap(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "acg";
  const char* seq_b = "tttacgttt"; 

  int match = 1;
  int mismatch = -1;
  int gap_open = -4;
  int gap_extend = -1;
  
  char gaps_ends_b = 0, no_gaps_in_a = 0, no_gaps_in_b = 0, freestartgap = 1;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               freestartgap, 0,
               no_gaps_in_a, no_gaps_in_b, 0, gaps_ends_b, 0, 0);
               
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_a, "------acg") == 0 && strcmp(result->result_b, "tttacgttt")== 0);
}

/* First sequence is aligned to the start of the second one because of the reduced cost of
 * mismatches and cost-free end gap */
void test_free_end_gap(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "acg";
  const char* seq_b = "tttacgttt"; 

  int match = 1;
  int mismatch = -1;
  int gap_open = -4;
  int gap_extend = -1;
  
  char gaps_ends_b = 0, no_gaps_in_a = 0, no_gaps_in_b = 0, freeendgap = 1;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               0, freeendgap,
               no_gaps_in_a, no_gaps_in_b, 0, gaps_ends_b, 0, 0);
               
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  TEST_CHECK(strcmp(result->result_a, "acg---------") == 0 && strcmp(result->result_b, "---tttacgttt")== 0);
}

void test(void)
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);
  
  const char* seq_a = "acg";
  const char* seq_b = "tttacgttt"; 

  int match = 1;
  int mismatch = -1;
  int gap_open = -4;
  int gap_extend = -1;
  
  char gaps_ends_b = 0, no_gaps_in_a = 0, no_gaps_in_b = 1, freeendgap = 1;
  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               1, freeendgap,
               no_gaps_in_a, no_gaps_in_b, 0, gaps_ends_b, 0, 0);
               
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  alignment_print_matrices(nw);
  TEST_CHECK(strcmp(result->result_a, "---acg---") == 0 && strcmp(result->result_b, "tttacgttt")== 0);
}


TEST_LIST = {
  { "no_gaps_in_longer", test_no_gaps_in_longer },
  { "no_gaps_equal_length", test_no_gaps_equal},
  { "gaps_only_at_ends_in_shorter", test_gaps_only_at_ends_in_shorter},
  { "gaps_only_at_ends_in_longer", test_gaps_only_at_ends_in_longer},
  { "no_gaps_smith_waterman", test_no_gaps_smith_waterman},
  { "free_start_gap", test_free_start_gap },
  { "free_end_gap", test_free_end_gap },
  { "test", test },
  { 0 }
};