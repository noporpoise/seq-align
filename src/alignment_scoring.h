/*
 alignment_scoring.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#ifndef ALIGNMENT_SCORING_HEADER_SEEN
#define ALIGNMENT_SCORING_HEADER_SEEN

#include <inttypes.h>
#include <stdbool.h>

typedef int score_t;

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

#define scoring_get_bit(arr,i) (((arr)[(i)/32] >> ((i)%32))&0x1)
#define scoring_set_bit(arr,i) ((arr)[(i)/32] |= (0x1<<((i)%32)))

#define get_wildcard_bit(scoring,c) scoring_get_bit((scoring)->wildcards,c)
#define set_wildcard_bit(scoring,c) scoring_set_bit((scoring)->wildcards,c)

#define get_swap_bit(scoring,a,b) scoring_get_bit((scoring)->swap_set[(size_t)a],b)
#define set_swap_bit(scoring,a,b) scoring_set_bit((scoring)->swap_set[(size_t)a],b)


void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool no_start_gap_penalty, bool no_end_gap_penalty,
                  bool no_gaps_in_a, bool no_gaps_in_b,
                  bool no_mismatches, bool case_sensitive);

/*
enum scoring_flags {
  NO_LEFT_GAP_PENALTY_A = 1, NO_RGHT_GAP_PENALTY_A = 2,
  NO_LEFT_GAP_A = 4, NO_GAPS_IN_A = 8, NO_RGHT_GAP_A = 16,
  NO_LEFT_GAP_PENALTY_B = 32, NO_RGHT_GAP_PENALTY_B = 64,
  NO_LEFT_GAP_B = 128, NO_GAPS_IN_B = 256, NO_RGHT_GAP_B = 512,
  NO_MISMATCHES = 1024, CASE_INSENSITIVE = 2048
};

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

#endif
