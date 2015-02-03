/*
 alignment_scoring.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#ifndef ALIGNMENT_SCORING_HEADER_SEEN
#define ALIGNMENT_SCORING_HEADER_SEEN

#include <inttypes.h>

typedef int score_t;

typedef struct
{
  int gap_open, gap_extend;

  // Needleman Wunsch only
  // Turn these on to turn off penalties for gaps at the start/end of alignment
  char no_start_gap_penalty, no_end_gap_penalty;

  // Turn at most one of these on at a time to prevent gaps/mismatches
  char no_gaps_in_a, no_gaps_in_b, no_mismatches;
  
  char gaps_ends_a, gaps_ends_b;

  // If swap_score not set, should we use match/mismatch values?
  char use_match_mismatch;
  int match, mismatch;

  char case_sensitive;

  // Array of characters that match to everything with the same penalty (i.e. 'N's)
  uint32_t wildcards[256/32], swap_set[256][256/32];
  score_t wildscores[256], swap_scores[256][256];
} scoring_t;

void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  char no_start_gap_penalty, char no_end_gap_penalty,
                  char no_gaps_in_a, char no_gaps_in_b, char gaps_ends_a, char gaps_ends_b,
                  char no_mismatches, char case_sensitive);

void scoring_add_wildcard(scoring_t* scoring, char c, int s);

void scoring_add_mutation(scoring_t* scoring, char a, char b, int score);

void scoring_print(const scoring_t* scoring);

void scoring_lookup(const scoring_t* scoring, char a, char b,
                    int* score, char* is_match);

// Some scoring systems
void scoring_system_PAM30(scoring_t *scoring);
void scoring_system_PAM70(scoring_t *scoring);
void scoring_system_BLOSUM80(scoring_t *scoring);
void scoring_system_BLOSUM62(scoring_t *scoring);
void scoring_system_DNA_hybridization(scoring_t *scoring);
void scoring_system_default(scoring_t *scoring); // DNA/RNA default
void scoring_system_lcs(scoring_t *scoring);

#endif
