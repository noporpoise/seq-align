/*
 tests/tests.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Feb 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h> // needed for rand()
#include <unistd.h>  // need for getpid() for getting setting rand number


#include "needleman_wunsch.h"
#include "smith_waterman.h"

//
// Tests
//
const char *suite_name = NULL;
char suite_pass = 0;
int suites_run = 0, suites_failed = 0, suites_empty = 0;
int tests_in_suite = 0, tests_run = 0, tests_failed = 0;

#define QUOTE(str) #str
#define ASSERT(x) {tests_run++; tests_in_suite++; if(!(x)) \
  { fprintf(stderr, "failed assert [%s:%i] %s\n", __FILE__, __LINE__, QUOTE(x)); \
    suite_pass = 0; tests_failed++; }}

void SUITE_START(const char *name)
{
  suite_pass = 1;
  suite_name = name;
  suites_run++;
  tests_in_suite = 0;
}

void SUITE_END()
{
  printf("Testing %s ", suite_name);
  size_t suite_i;
  for(suite_i = strlen(suite_name); suite_i < 80-8-5; suite_i++) printf(".");
  printf("%s\n", suite_pass ? " pass" : " fail");
  if(!suite_pass) suites_failed++;
  if(!tests_in_suite) suites_empty++;
}

//
// Useful MACROs

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0);

#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MIN(x,y) ((x) <= (y) ? (x) : (y))

//

/* No gap is expected in the longer sequence */
void nw_test_no_gaps_in_longer()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);

  const char* seq_a = "aaaaacg";
  const char* seq_b = "acgt";

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = true, no_gaps_in_b = false;
  bool no_mismatches = false, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, aln);

  // ASSERT(strcmp(aln->result_a, "aaaaacg") == 0 &&
  //        strcmp(aln->result_b, "acgt---") == 0);

  ASSERT(strcmp(aln->result_a, "aaaaacg-") == 0 &&
         strcmp(aln->result_b, "a----cgt") == 0);

  alignment_free(aln);
  needleman_wunsch_free(nw);
}

/* First sequence is aligned to the corresponding (equal) substring of the second
 * sequence because both gaps at start and at end are free */
void nw_test_free_gaps_at_ends()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  const char* seq_a = "acg";
  const char* seq_b = "tttacgttt";

  int match = 1;
  int mismatch = -1;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = true, no_end_gap_penalty = true;
  bool no_gaps_in_a = false, no_gaps_in_b = false;
  bool no_mismatches = false, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);
  ASSERT(strcmp(result->result_a, "---acg---") == 0 &&
         strcmp(result->result_b, "tttacgttt") == 0);

  alignment_free(result);
  needleman_wunsch_free(nw);
}

void nw_test_no_mismatches()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *result = alignment_create(256);

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = false, no_gaps_in_b = false;
  bool no_mismatches = true, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  needleman_wunsch_align("atc", "ac", &scoring, nw, result);
  ASSERT(strcmp(result->result_a, "atc") == 0 &&
         strcmp(result->result_b, "a-c") == 0);

  needleman_wunsch_align("cgatcga", "catcctcga", &scoring, nw, result);
  ASSERT(strcmp(result->result_a, "cgatc---ga") == 0 &&
         strcmp(result->result_b, "c-atcctcga") == 0);

  alignment_free(result);
  needleman_wunsch_free(nw);
}

// size is buffer size including NUL byte
static char* make_rand_seq(char *seq, size_t size)
{
  if(size == 0) { return seq; }
  size_t i, len = rand() % (size-1);
  const char bases[] = "acgt";
  for(i = 0; i < len; i++) seq[i] = bases[rand() & 3];
  seq[len] = '\0';
  return seq;
}

void nw_test_no_mismatches_rand()
{
  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = false, no_gaps_in_b = false;
  bool no_mismatches = true, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  char seqa[100], seqb[100];
  size_t i;

  // Run 50 random alignments
  for(i = 0; i < 50; i++)
  {
    make_rand_seq(seqa, sizeof(seqa));
    make_rand_seq(seqb, sizeof(seqb));
    needleman_wunsch_align(seqa, seqb, &scoring, nw, aln);
    // Check no mismatches
    char *a = aln->result_a, *b = aln->result_b;
    while(1) {
      ASSERT(*a == '-' || *b == '-' || *a == *b);
        // printf("Seq: '%s' '%s'\n", aln->result_a, aln->result_b);
        // exit(EXIT_FAILURE);
      if(!*a && !*b) break;
      a++; b++;
    }
  }

  alignment_free(aln);
  needleman_wunsch_free(nw);
}


void test_nw()
{
  SUITE_START("Needleman-Wunsch");

  nw_test_no_gaps_in_longer();
  nw_test_free_gaps_at_ends();
  nw_test_no_mismatches();
  nw_test_no_mismatches_rand();

  SUITE_END();
}

void sw_test_no_gaps_smith_waterman()
{
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *result = alignment_create(256);

  const char* seq_a = "gacag";
  const char* seq_b = "tgaagt";

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool no_start_gap_penalty = false, no_end_gap_penalty = false;
  bool no_gaps_in_a = true, no_gaps_in_b = true;
  bool no_mismatches = false, case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b,
               no_mismatches, case_sensitive);

  smith_waterman_align(seq_a, seq_b, &scoring, sw);

  smith_waterman_fetch(sw, result);
  ASSERT(strcmp(result->result_a, "ga") == 0 &&
         strcmp(result->result_b, "ga") == 0);

  smith_waterman_fetch(sw, result);
  ASSERT(strcmp(result->result_a, "ag") == 0 &&
         strcmp(result->result_b, "ag") == 0);

  alignment_free(result);
  smith_waterman_free(sw);
}


void test_sw()
{
  SUITE_START("Smith-Waterman");

  sw_test_no_gaps_smith_waterman();

  SUITE_END();
}

int main(int argc, char **argv)
{
  if(argc != 1)
  {
    printf("  Unused args '%s..'\n", argv[1]);
    printf("Usage: ./bit_array_test\n");
    exit(EXIT_FAILURE);
  }

  // Initialise random number generator
  srand((unsigned int)time(NULL) + getpid());

  printf("  Test seq-align C library:\n\n");

  // Test suites go here
  test_nw();
  test_sw();

  printf("\n");
  printf(" %i / %i suites failed\n", suites_failed, suites_run);
  printf(" %i / %i suites empty\n", suites_empty, suites_run);
  printf(" %i / %i tests failed\n", tests_failed, tests_run);

  printf("\n THE END.\n");

  return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
