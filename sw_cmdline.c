/*
 sw_cmdline.c
 project: SmithWaterman
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/smithwaterman/
 Copyright (C) 18-Dec-2011

 see README

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define MAX(a,b)   ((a) >= (b) ? (a) : (b))
#define MIN(a,b)   ((a) <= (b) ? (a) : (b))

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h> // required for va_list

// my utility functions
#include "string_buffer.h"
#include "seq_file.h"
#include "utility_lib.h"

// Alignment scoring and loading
#include "alignment_scoring_load.h"
#include "alignment_cmdline.h"

#include "smith_waterman.h"

// For this run
char* cmd;
char print_colour = 0, print_seq = 0, print_fasta = 0, print_pretty = 0;
char interactive = 0;

score_t min_score = 0;
char min_score_set = 0;

unsigned long max_hits_per_alignment = 0;
char max_hits_per_alignment_set = 0;

unsigned int print_context = 0;
SCORING_SYSTEM* scoring;

unsigned long alignment_index = 0;

void set_default_scoring()
{
  scoring = scoring_system_default();

  // Change slightly
  scoring->match = 2;
  scoring->mismatch = -2;
  scoring->gap_open = -2;
  scoring->gap_open = -1;
}

void print_usage(char* err_fmt, ...)
{
  if(err_fmt != NULL)
  {
    va_list argptr;
    va_start(argptr, err_fmt);
    
    StrBuf *error = strbuf_init(200);
    strbuf_append_str(error, "SmithWaterman Error: ");
    strbuf_vsprintf(error, strbuf_len(error), err_fmt, argptr);
    strbuf_chomp(error);

    va_end(argptr);

    fprintf(stderr, "%s\n", error->buff);
    fprintf(stderr, "Use -h option to print help\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  if(scoring != NULL)
  {
    scoring_free(scoring);
  }

  // Get and print defaults
  set_default_scoring();

  fprintf(stderr, "usage: %s [OPTIONS] [seq1 seq2]\n", cmd);

  fprintf(stderr,
"  Smith-Waterman optimal local alignment (maximises score).  \n"
"  Takes a pair of sequences on the command line, or can read from a\n"
"  file and from sequence piped in.  Can read gzip files and FASTA.\n"
"\n");

  fprintf(stderr,
          "  OPTIONS:\n"
"    --file <file>        Sequence file reading with gzip support - read two\n"
"                         sequences at a time and align them\n"
"    --files <f1> <f2>    Read one sequence from each file to align at one time\n"
"    --stdin              Read from STDIN (same as '--file -')\n"
"\n"
"    --case_sensitive     Use case sensitive character comparison [default: off]\n"
"\n");

  fprintf(stderr,
"    --match <score>      [default: %i]\n"
"    --mismatch <score>   [default: %i]\n"
"    --gapopen <score>    [default: %i]\n"
"    --gapextend <score>  [default: %i]\n"
"\n"
"    --nogaps             No gaps allowed in the alignment\n"
"    --nomismatches       No mismatches allowed. \n"
"      When used together, prints longest common substrings in order of length\n"
"\n"
"    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>\n"
"    --substitution_matrix <file>  see details for formatting\n"
"    --substitution_pairs <file>   see details for formatting\n"
"\n"
"    --wildcard <w> <s>   Character <w> matches all characters with score <s>\n"
"\n",
          scoring->match, scoring->mismatch,
          scoring->gap_open, scoring->gap_extend);

  fprintf(stderr,
"    --minscore <score>   Minimum required score\n"
"                         [default: match * MAX(0.2 * length, 2)]\n"
"    --maxhits <hits>     Maximum number of results per alignment\n"
"                         [default: no limit]\n"
"    --context <n>        Print <n> bases of context\n"
"    --printfasta         Print fasta header lines\n"
"    --printseq           Print sequences before local alignments\n"
"    --pretty             Print with a descriptor line\n"
"    --colour             Print with colour\n"
"\n"
" DETAILS:\n"
"  * For help choosing scoring, see the README file. \n"
"  * Gap (of length N) penalty is: (open+N*extend)\n"
//"  * To do alignment without affine gap, set '--gapopen 0'.\n"
"  * Scoring files should be matrices, with entries separated by a single\n"
"    character or whitespace. See files in the 'scores' directory for examples.\n"
"\n"
"  turner.isaac@gmail.com  (compiled: %s %s)\n", __DATE__, __TIME__);

  exit(EXIT_FAILURE);
}

// Print one line of an alignment
void print_alignment_part(const char* seq1, const char* seq2,
                          size_t pos, size_t len,
                          const char* context_str,
                          size_t spaces_left, size_t spaces_right,
                          size_t context_left, size_t context_right)
{
  printf("  ");

  unsigned int i;
  for(i = 0; i < spaces_left; i++)
    printf(" ");

  if(context_left > 0)
  {
    if(print_colour)
    {
      printf("%s", align_col_context);
    }

    printf("%.*s", (int)context_left, context_str+pos-context_left);

    if(print_colour)
    {
      printf("%s", align_col_stop);
    }
  }

  if(print_colour)
  {
    alignment_colour_print_against(seq1, seq2, scoring->case_sensitive);
  }
  else
  {
    printf("%s", seq1);
  }
  
  if(context_right > 0)
  {
    if(print_colour)
    {
      printf("%s", align_col_context);
    }

    printf("%.*s", (int)context_right, context_str+pos+len);

    if(print_colour)
    {
      printf("%s", align_col_stop);
    }
  }

  for(i = 0; i < spaces_right; i++)
    printf(" ");

  printf("  [pos: %li; len: %lu]\n", pos, len);
}

char get_next_hit()
{
  if(!interactive)
    return 1;

  int r = 0;

  char response = 0;
  char next_hit = 0;

  while(!response)
  {
    printf("next [h]it or [a]lignment: ");
    fflush(stdout);

    while((r = getc(stdin)) != -1 && r != '\n' && r != '\r')
    {
      if(r == 'h' || r == 'H')
      {
        next_hit = 1;
        response = 1;
      }
      else if(r == 'a' || r == 'A')
      {
        next_hit = 0;
        response = 1;
      }
    }

    if(r == -1)
    {
      // We're done here
      putc('\n', stdout);
      exit(EXIT_SUCCESS);
    }
  }

  return next_hit;
}

// Align two sequences against each other to find local alignments between them
void align(const char *seq_a, const char *seq_b,
           const char *seq_a_name, const char *seq_b_name)
{
  // Check both arguments have length > 0
  if(seq_a[0] == '\0' || seq_b[0] == '\0')
  {
    fprintf(stderr, "SmithWaterman Error: Sequences must have length > 0\n");
    fflush(stderr);

    if(print_fasta && seq_a_name != NULL && seq_b_name != NULL)
    {
      fprintf(stderr, "%s\n%s\n",seq_a_name,seq_b_name);
    }

    fflush(stderr);

    return;
  }

  SW_COMPUTATION* smithwaterman = smith_waterman_align(seq_a, seq_b, scoring);

  size_t len_a = smith_waterman_seq_a_strlen(smithwaterman);
  size_t len_b = smith_waterman_seq_b_strlen(smithwaterman);

  printf("== Alignment %lu lengths (%lu, %lu):\n", alignment_index,
         (unsigned long)len_a, (unsigned long)len_b);

  // seqA
  if(print_fasta && seq_a_name != NULL)
  {
    printf("%s\n", seq_a_name);
  }
  
  if(print_seq)
  {
    printf("%s\n", seq_a);
  }
  
  // seqB
  if(print_fasta && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_seq)
  {
    printf("%s\n", seq_b);
  }

  printf("\n");

  if(!min_score_set)
  {
    // If min_score hasn't been set, set a limit based on the lengths of seqs
    // or zero if we're running interactively
    min_score = interactive ? 0 : scoring->match * MAX(0.2 * MIN(len_a, len_b), 2);

    #ifdef DEBUG
    printf("min_score: %i\n", min_score);
    #endif
  }

  fflush(stdout);

  // Allocate memory for storing result
  SW_LOCAL_ALIGNMENT* alignment = smith_waterman_create_hit();

  unsigned long hit_index = 0;

  // For print context
  size_t context_left = 0, context_right = 0;
  size_t left_spaces_a = 0, left_spaces_b = 0;
  size_t right_spaces_a = 0, right_spaces_b = 0;


  while(get_next_hit() &&
        smith_waterman_get_hit(smithwaterman, alignment) &&
        alignment->score >= min_score &&
        (!max_hits_per_alignment_set || hit_index < max_hits_per_alignment))
  {
    printf("hit %lu.%lu score: %i\n", alignment_index, hit_index++,
           alignment->score);
    
    if(print_context)
    {
      // Calculate number of characters of context to print either side
      context_left = MAX(alignment->pos_a, alignment->pos_b);
      context_left = MIN(context_left, print_context);

      size_t rem_a = len_a - (alignment->pos_a + alignment->len_a);
      size_t rem_b = len_b - (alignment->pos_b + alignment->len_b);

      context_right = MAX(rem_a, rem_b);
      context_right = MIN(context_right, print_context);
    
      left_spaces_a = (context_left > alignment->pos_a)
                      ? context_left - alignment->pos_a : 0;

      left_spaces_b = (context_left > alignment->pos_b)
                      ? context_left - alignment->pos_b : 0;
    
      right_spaces_a = (context_right > rem_a) ? context_right - rem_a : 0;
      right_spaces_b = (context_right > rem_b) ? context_right - rem_b : 0;
    }

    #ifdef DEBUG
    printf("context left = %lu; right = %lu spacing: [%lu,%lu] [%lu,%lu]\n",
           context_left, context_right,
           left_spaces_a, right_spaces_a,
           left_spaces_b, right_spaces_b);
    #endif

    // seq a
    print_alignment_part(alignment->result_a, alignment->result_b,
                         alignment->pos_a, alignment->len_a,
                         seq_a,
                         left_spaces_a, right_spaces_a,
                         context_left-left_spaces_a,
                         context_right-right_spaces_a);

    if(print_pretty)
    {
      printf("  ");

      size_t max_left_spaces = MAX(left_spaces_a, left_spaces_b);
      size_t max_right_spaces = MAX(right_spaces_a, right_spaces_b);

      size_t spacer;

      // Print spaces for lefthand spacing
      for(spacer = 0; spacer < max_left_spaces; spacer++)
      {
        printf(" ");
      }

      // Print dots for lefthand context sequence
      for(spacer = 0; spacer < context_left-max_left_spaces; spacer++)
      {
        printf(".");
      }
      
      alignment_print_spacer(alignment->result_a, alignment->result_b, scoring);

      // Print dots for righthand context sequence
      for(spacer = 0; spacer < context_right-max_right_spaces; spacer++)
      {
        printf(".");
      }

      // Print spaces for righthand spacing
      for(spacer = 0; spacer < max_right_spaces; spacer++)
      {
        printf(" ");
      }

      printf("\n");
    }

    // seq b
    print_alignment_part(alignment->result_b, alignment->result_a,
                         alignment->pos_b, alignment->len_b,
                         seq_b,
                         left_spaces_b, right_spaces_b,
                         context_left-left_spaces_b,
                         context_right-right_spaces_b);

    printf("\n");

    // Flush output here
    fflush(stdout);
  }

  printf("==\n");
  fflush(stdout);

  // Free memory used for holding results
  smith_waterman_free_hit(alignment);

  // Free alignment between the two sequences
  smith_waterman_free(smithwaterman);

  // Increment sequence alignment counter
  alignment_index++;
}

void align_pair_from_file(StrBuf* seq_a, StrBuf *seq_b,
                          const char *seq_a_name, const char *seq_b_name)
{
  align(seq_a->buff, seq_b->buff, seq_a_name, seq_b_name);
}

int main(int argc, char* argv[])
{
  cmd = argv[0];

  #ifdef DEBUG
  printf("DEBUG: on\n");
  #endif
  
  if(argc == 1)
  {
    print_usage(NULL);
  }
  
  //
  // Command line arguments handled here
  //

  char *seq1 = NULL, *seq2 = NULL;

  cmdline_init();

  scoring = NULL;
  char case_sensitive = 0;

  // First run through arguments to set up case_sensitive and scoring system

  // case sensitive needs to be dealt with first
  // (it is used to construct hash table for swap_table)
  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(strcasecmp(argv[argi], "--help") == 0 ||
       strcasecmp(argv[argi], "-help") == 0 ||
       strcasecmp(argv[argi], "-h") == 0)
    {
      print_usage(NULL);
    }
    else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
    {
      case_sensitive = 1;
    }
    else if(strcasecmp(argv[argi], "--scoring") == 0)
    {
      if(scoring != NULL)
      {
        print_usage("More than one scoring system specified - not permitted");
      }
    
      if(strcasecmp(argv[argi+1], "PAM30") == 0)
      {
        scoring = scoring_system_PAM30();
      }
      else if(strcasecmp(argv[argi+1], "PAM70") == 0)
      {
        scoring = scoring_system_PAM70();
      }
      else if(strcasecmp(argv[argi+1], "BLOSUM80") == 0)
      {
        scoring = scoring_system_BLOSUM80();
      }
      else if(strcasecmp(argv[argi+1], "BLOSUM62") == 0)
      {
        scoring = scoring_system_BLOSUM62();
      }
      else if(strcasecmp(argv[argi+1], "DNA_HYBRIDIZATION") == 0)
      {
        scoring = scoring_system_DNA_hybridization();
      }
      else {
        print_usage("Unknown --scoring choice, not one of "
                    "PAM30|PAM70|BLOSUM80|BLOSUM62");
      }

      argi++; // took an argument
    }
  }

  // Set up default scoring now
  if(scoring == NULL)
  {
    set_default_scoring();
  }

  scoring->case_sensitive = case_sensitive;
  // Scoring is now initiated - may tweak later

  // Keep track of what is set
  char substitutions_set = 0;
  char match_set = 0;
  char mismatch_set = 0;

  for(argi = 1; argi < argc; argi++)
  {
    if(argv[argi][0] == '-')
    {
      // strcasecmp does case insensitive comparison
      if(strcasecmp(argv[argi], "--nogaps") == 0)
      {
        scoring->no_gaps_in_a = 1;
        scoring->no_gaps_in_b = 1;
      }
      else if(strcasecmp(argv[argi], "--nogapsin1") == 0)
      {
        scoring->no_gaps_in_a = 1;
      }
      else if(strcasecmp(argv[argi], "--nogapsin2") == 0)
      {
        scoring->no_gaps_in_b = 1;
      }
      else if(strcasecmp(argv[argi], "--nomismatches") == 0)
      {
        scoring->no_mismatches = 1;
      }
      else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
      {
        // Already dealt with
        //case_sensitive = 1;
      }
      else if(strcasecmp(argv[argi], "--printseq") == 0)
      {
        print_seq = 1;
      }
      else if(strcasecmp(argv[argi], "--printfasta") == 0)
      {
        print_fasta = 1;
      }
      else if(strcasecmp(argv[argi], "--pretty") == 0)
      {
        print_pretty = 1;
      }
      else if(strcasecmp(argv[argi], "--colour") == 0)
      {
        print_colour = 1;
      }
      else if(strcasecmp(argv[argi], "--stdin") == 0)
      {
        // Similar to --file argument below
        cmdline_add_files("-", NULL);
        interactive = 1;
      }
      else if(argi == argc-1)
      {
        // All the remaining options take an extra argument
        print_usage("Unknown argument without parameter %s",argv[argi]);
      }
      else if(strcasecmp(argv[argi], "--scoring") == 0)
      {
        // This handled above
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_matrix") == 0)
      {
        gzFile* sub_matrix_file = gzopen(argv[argi+1], "r");
        // gzbuffer(sub_matrix_file, 16384); // doesn't seem to work

        align_scoring_load_matrix(sub_matrix_file, argv[argi+1],
                                  scoring, case_sensitive);

        //gzclose_r(sub_matrix_file); // doesn't seem to work
        gzclose(sub_matrix_file);
        substitutions_set = 1;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_pairs") == 0)
      {
        gzFile* sub_pairs_file = gzopen(argv[argi+1], "r");
        //gzbuffer(sub_pairs_file, 16384); // doesn't seem to work
        
        align_scoring_load_pairwise(sub_pairs_file, argv[argi+1],
                                    scoring, case_sensitive);
        
        //gzclose_r(sub_pairs_file); // doesn't seem to work
        gzclose(sub_pairs_file);
        substitutions_set = 1;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--minscore") == 0)
      {
        if(!parse_entire_uint(argv[argi+1], &min_score))
        {
          print_usage("Invalid --minscore <score> argument (must be a +ve int)");
        }
        min_score_set = 1;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--maxhits") == 0)
      {
        if(!parse_entire_ulong(argv[argi+1], &max_hits_per_alignment))
        {
          print_usage("Invalid --maxhits <hits> argument (must be a +ve int)");
        }
        max_hits_per_alignment_set = 1;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--context") == 0)
      {
        if(!parse_entire_uint(argv[argi+1], &print_context))
        {
          print_usage("Invalid --context <c> argument (must be >= 0)");
        }

        argi++;
      }
      else if(strcasecmp(argv[argi], "--match") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->match))
        {
          print_usage("Invalid --match argument ('%s') must be an int",
                      argv[argi+1]);
        }

        match_set = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--mismatch") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->mismatch))
        {
          print_usage("Invalid --mismatch <penalty> argument (must be an int)");
        }

        mismatch_set = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_open))
        {
          print_usage("Invalid --gapopen <penalty> argument (must be an int)");
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_extend))
        {
          print_usage("Invalid --gapextend <penalty> argument (must be an int)");
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--file") == 0)
      {
        cmdline_add_files(argv[argi+1], NULL);
        argi++; // took an argument
      }
      // Remaining options take two arguments but check themselves
      else if(strcasecmp(argv[argi], "--files") == 0)
      {
        if(argi >= argc-2)
        {
          print_usage("--files option takes 2 arguments");
        }
        else if(strcmp(argv[argi+1], "-") == 0 &&
                strcmp(argv[argi+2], "-") == 0)
        {
          // Read both from stdin
          cmdline_add_files(argv[argi+1], NULL);
        }
        else
        {
          cmdline_add_files(argv[argi+1], argv[argi+2]);
        }

        argi += 2; // took two arguments
      }
      else if(strcasecmp(argv[argi], "--wildcard") == 0)
      {
        int wildscore = 0;

        if(argi == argc-2 || strlen(argv[argi+1]) != 1 ||
           !parse_entire_int(argv[argi+2], &wildscore))
        {
          print_usage("--wildcard <w> <s> takes a single character and a number");
        }

        scoring_add_wildcard(scoring, argv[argi+1][0], wildscore);

        argi += 2; // took two arguments
      }
      else
      {
        // Error - unknown option
        print_usage("Unknown argument '%s'",argv[argi]);
      }
    }
    else
    {
      if(argc - argi != 2)
      {
        print_usage("Unknown argument(s) '%s'", argv[argi]);
      }
      break;
    }
  }

  if((match_set && !mismatch_set && !scoring->no_mismatches) ||
     (!match_set && mismatch_set))
  {
    print_usage("--match --mismatch must both be set or neither set");
  }
  else if(substitutions_set && !match_set)
  {
    // if substitution table set and not match/mismatch
    scoring->use_match_mismatch = 0;
  }

  // Check for extra unused arguments
  // and set seq1 and seq2 if they have been passed
  if(argi < argc)
  {
    seq1 = argv[argi];
    seq2 = argv[argi+1];
  }

  int file_list_length = cmdline_get_num_of_file_pairs();

  if(seq1 == NULL && file_list_length == 0)
  {
    print_usage("No input specified");
  }

  // ALIGN!
  if(seq1 != NULL)
  {
    // Align seq1 and seq2
    align(seq1, seq2, NULL, NULL);
  }

  int i;
  for(i = 0; i < file_list_length; i++)
  {
    // Align from files
    align_from_file(cmdline_get_file1(i), cmdline_get_file2(i),
                    &align_pair_from_file);
  }

  cmdline_finish();

  scoring_free(scoring);

  return EXIT_SUCCESS;
}
