/*
 alignment_scoring.h
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

#ifndef ALIGNMENT_SCORING_HEADER_SEEN
#define ALIGNMENT_SCORING_HEADER_SEEN

#include <inttypes.h>
#include <stdbool.h>
#include <limits.h> // INT_MIN

#include "bit_array/bit_macros.h"

typedef int score_t;
#define SCORE_MIN INT_MIN

typedef struct
{
  int gap_open, gap_extend;

  // Needleman Wunsch only
  // Turn these on to turn off penalties for gaps at the start/end of alignment
  bool no_start_gap_penalty, no_end_gap_penalty;

  // Turn at most one of these on at a time to prevent gaps/mismatches
  bool no_gaps_in_a, no_gaps_in_b, no_mismatches;

  // If swap_score not set, should we use match/mismatch values?
  bool use_match_mismatch;
  int match, mismatch;

  bool case_sensitive;

  // Array of characters that match to everything with the same penalty (i.e. 'N's)
  uint32_t wildcards[256/32], swap_set[256][256/32];
  score_t wildscores[256], swap_scores[256][256];
} scoring_t;

#define get_wildcard_bit(scoring,c) bitset_get((scoring)->wildcards,c)
#define set_wildcard_bit(scoring,c) bitset_set((scoring)->wildcards,c)

#define get_swap_bit(scoring,a,b) bitset_get((scoring)->swap_set[(size_t)a],b)
#define set_swap_bit(scoring,a,b) bitset_set((scoring)->swap_set[(size_t)a],b)

#ifdef __cplusplus
extern "C" {
#endif

void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool no_start_gap_penalty, bool no_end_gap_penalty,
                  bool no_gaps_in_a, bool no_gaps_in_b,
                  bool no_mismatches, bool case_sensitive);

/*

enum scoring_flags
{
  FREE_LEFT_GAP_A         = 0x1, FREE_LEFT_GAP_B         = 0x2,
  FREE_LEFT_GAP_SHORTEST  = 0x4, FREE_LEFT_GAP_LONGEST   = 0x8,

  FREE_RGHT_GAP_A         = 0x10, FREE_RGHT_GAP_B        = 0x20,
  FREE_RGHT_GAP_SHORTEST  = 0x40, FREE_RGHT_GAP_LONGEST  = 0x80,

  NO_LEFT_GAP_A           = 0x100, NO_LEFT_GAP_B         = 0x200,
  NO_LEFT_GAP_SHORTEST    = 0x400, NO_LEFT_GAP_LONGEST   = 0x800,

  NO_GAPS_IN_A            = 0x1000, NO_GAPS_IN_B         = 0x2000,
  NO_GAPS_IN_SHORTEST     = 0x4000, NO_GAPS_IN_LONGEST   = 0x8000,

  NO_RGHT_GAP_A           = 0x10000, NO_RGHT_GAP_B       = 0x10000,
  NO_RGHT_GAP_SHORTEST    = 0x40000, NO_RGHT_GAP_LONGEST = 0x80000,

  NO_MISMATCHES = 0x100000, CASE_INSENSITIVE = 0x200000
};

typedef struct
{
  int match, mismatch;
  int gap_open, gap_extend;
  int flags;

  // Array of characters that match to everything with the same penalty (i.e. 'N's)
  uint32_t wildcards[256/32], swap_set[256][256/32];
  score_t wildscores[256], swap_scores[256][256];
} scoring_t;

void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool no_start_gap_penalty, bool no_end_gap_penalty,
                  bool no_gaps_in_a, bool no_gaps_in_b,
                  bool no_mismatches, bool case_sensitive)
{
  scoring_init2(scoring, match, mismatch, gap_open, gap_extend,
                (no_start_gap_penalty ? FREE_LEFT_GAP_A | FREE_LEFT_GAP_B : 0) |
                (no_end_gap_penalty   ? FREE_RGHT_GAP_A | FREE_RGHT_GAP_B : 0) |
                (no_gaps_in_a ? NO_GAPS_IN_A : 0) |
                (no_gaps_in_b ? NO_GAPS_IN_B : 0) |
                (no_mismatches ? NO_MISMATCHES : 0) |
                (case_sensitive ? 0 : CASE_INSENSITIVE))
}

void scoring_init2(scoring_t* scoring,
                   int match, int mismatch,
                   int gap_open, int gap_extend,
                   int flags);
*/

#define scoring_is_wildcard(scoring,c) (get_wildcard_bit(scoring,c))

void scoring_add_wildcard(scoring_t* scoring, char c, int s);

void scoring_add_mutation(scoring_t* scoring, char a, char b, int score);

void scoring_print(const scoring_t* scoring);

void scoring_lookup(const scoring_t* scoring, char a, char b,
                    int *score, bool *is_match);

// Some scoring systems
void scoring_system_PAM30(scoring_t *scoring);
void scoring_system_PAM70(scoring_t *scoring);
void scoring_system_BLOSUM80(scoring_t *scoring);
void scoring_system_BLOSUM62(scoring_t *scoring);
void scoring_system_DNA_hybridization(scoring_t *scoring);
void scoring_system_default(scoring_t *scoring); // DNA/RNA default

#ifdef __cplusplus
}
#endif

#endif
