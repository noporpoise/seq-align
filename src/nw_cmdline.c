/*
 nw_cmdline.c
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower
#include <stdarg.h> // required for va_list

// my utility functions
#include "string_buffer.h"
#include "seq_file.h"

// Alignment scoring and loading
#include "alignment_scoring_load.h"
#include "alignment_cmdline.h"

#include "needleman_wunsch.h"

// For this run
char* cmd;
char print_colour = 0, print_pretty = 0, print_scores = 0,
     print_fasta = 0, print_zam = 0;

scoring_t scoring;

// Alignment results stored here
nw_aligner_t *nw;
alignment_t *result;

void set_default_scoring()
{
  scoring_system_default(&scoring);
}

static void print_usage(const char* errfmt, ...) __attribute__((noreturn));
static void print_usage(const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fprintf(stderr, "\n");
  }

  // Get and print defaults
  set_default_scoring();

  fprintf(stderr, "usage: %s [OPTIONS] [seq1 seq2]\n", cmd);

  fprintf(stderr,
"  Needleman-Wunsch optimal global alignment (maximises score). Takes a pair \n"
"  of sequences on the command line, reads from a file and from sequence \n"
"  piped in.  Can read gzip files and those in FASTA, FASTQ or plain format.\n"
"\n"
"  OPTIONS:\n"
"    --file <file>        Sequence file reading with gzip support\n"
"    --files <f1> <f2>    Read one sequence from each file at a time to align\n"
"    --stdin              Read from STDIN (same as '--file -')\n"
"\n"
"    --case_sensitive     Case sensitive character comparison\n"
"\n");

  fprintf(stderr,
"    --match <score>      [default: %i]\n"
"    --mismatch <score>   [default: %i]\n"
"    --gapopen <score>    [default: %i]\n"
"    --gapextend <score>  [default: %i]\n"
"\n"
"    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>\n"
"    --substitution_matrix <file>  see details for formatting\n"
"    --substitution_pairs <file>   see details for formatting\n"
"\n"
"    --wildcard <w> <s>   Character <w> matches all characters with score <s>\n",
          scoring.match, scoring.mismatch,
          scoring.gap_open, scoring.gap_extend);

  fprintf(stderr,
"\n"
"    --freestartgap       No penalty for gap at start of alignment\n"
"    --freeendgap         No penalty for gap at end of alignment\n"
"\n"
"    --printscores        Print optimal alignment scores\n"
"    --printfasta         Print fasta header lines\n"
"    --pretty             Print with a descriptor line\n"
"    --colour             Print with colour\n"
"    --zam                A funky type of output\n"
"\n"
"  EXPERIMENTAL (and buggy):\n"
"    --nogapsin1          No gaps allowed in the first sequence\n"
"    --nogapsin2          No gaps allowed in the second sequence\n"
"    --nogaps             No gaps allowed in either sequence\n"
"    --nomismatches       No mismatches allowed: cannot be used with --nogaps..\n"
"\n"
" DETAILS:\n"
"  * For help choosing scoring, see the README file. \n"
"  * Gap (of length N) penalty is: (open+N*extend)\n"
//"  * To do alignment without affine gap, set '--gapopen 0'.\n"
"  * Scoring files should be matrices, with entries separated by a single \n"
"    character or whitespace.  See files in the 'scores' directory for examples.\n"
"\n"
"  turner.isaac@gmail.com  (compiled: %s %s)\n", __DATE__, __TIME__);

  exit(EXIT_FAILURE);
}

static void align_zam(const char *seq_a, const char *seq_b)
{
  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);

  // Swap '-' for '_'
  int i;
  for(i = 0; result->result_a[i] != '\0'; i++)
  {
    if(result->result_a[i] == '-')
    {
      result->result_a[i] = '_';
    }

    if(result->result_b[i] == '-')
    {
      result->result_b[i] = '_';
    }
  }

  int num_of_mismatches = 0;
  int num_of_indels = 0;

  // Print branch 1
  printf("Br1:%s\n", result->result_a);

  // Print spacer
  printf("    ");

  for(i = 0; result->result_a[i] != '\0'; i++)
  {
    if(result->result_a[i] == '_' || result->result_b[i] == '_')
    {
      printf(" ");
      num_of_indels++;
    }
    else if((scoring.case_sensitive && result->result_a[i] != result->result_b[i]) ||
            tolower(result->result_a[i]) != tolower(result->result_b[i]))
    {
      printf("*");
      num_of_mismatches++;
    }
    else
    {
      printf("|");
    }
  }

  printf("\n");

  // Print branch 2
  printf("Br2:%s\n", result->result_b);

  // print mismatch indel numbers
  printf("%i %i\n\n", num_of_mismatches, num_of_indels);
}

static void align(const char *seq_a, const char *seq_b,
                  const char *seq_a_name, const char *seq_b_name)
{
  if(print_zam)
  {
    align_zam(seq_a, seq_b);
    fflush(stdout);
    return;
  }

  needleman_wunsch_align(seq_a, seq_b, &scoring, nw, result);

  if(print_fasta && seq_a_name != NULL)
  {
    printf("%s\n", seq_a_name);
  }

  if(print_fasta && print_pretty && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_colour)
  {
    // Print alignment line 1
    alignment_colour_print_against(result->result_a, result->result_b,
                                   scoring.case_sensitive);
  }
  else
  {
    printf("%s", result->result_a);
  }
  printf("\n");
  
  if(print_pretty)
  {
    // Print spacer
    alignment_print_spacer(result->result_a, result->result_b, &scoring);

    printf("\n");
  }
  else if(print_fasta && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_colour)
  {
    // Print alignment line 2
    alignment_colour_print_against(result->result_b, result->result_a,
                                   scoring.case_sensitive);
  }
  else
  {
    printf("%s", result->result_b);
  }
  printf("\n");

  if(print_scores)
  {
    printf("score: %i\n", result->score);
  }
  
  printf("\n");

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

  char case_sensitive = 0;
  char scoring_set = 0;
  
  // First run through arguments to set up case_sensitive and scoring system

  // case sensitive needs to be dealt with first
  // (is is used to construct hash table for swap_scores)
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
      if(scoring_set)
      {
        print_usage("More than one scoring system specified - not permitted");
      }
    
      if(strcasecmp(argv[argi+1], "PAM30") == 0)
      {
        scoring_system_PAM30(&scoring);
      }
      else if(strcasecmp(argv[argi+1], "PAM70") == 0)
      {
        scoring_system_PAM70(&scoring);
      }
      else if(strcasecmp(argv[argi+1], "BLOSUM80") == 0)
      {
        scoring_system_BLOSUM80(&scoring);
      }
      else if(strcasecmp(argv[argi+1], "BLOSUM62") == 0)
      {
        scoring_system_BLOSUM62(&scoring);
      }
      else if(strcasecmp(argv[argi+1], "DNA_HYBRIDIZATION") == 0)
      {
        scoring_system_DNA_hybridization(&scoring);
      }
      else {
        print_usage("Unknown --scoring choice, not one of "
                    "PAM30|PAM70|BLOSUM80|BLOSUM62");
      }

      scoring_set = 1;
      argi++; // took an argument
    }
  }

  // Set up default scoring now
  if(!scoring_set)
  {
    set_default_scoring();
  }

  scoring.case_sensitive = case_sensitive;
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
      if(strcasecmp(argv[argi], "--freestartgap") == 0)
      {
        scoring.no_start_gap_penalty = 1;
      }
      else if(strcasecmp(argv[argi], "--freeendgap") == 0)
      {
        scoring.no_end_gap_penalty = 1;
      }
      else if(strcasecmp(argv[argi], "--nogaps") == 0)
      {
        scoring.no_gaps_in_a = 1;
        scoring.no_gaps_in_b = 1;
      }
      else if(strcasecmp(argv[argi], "--nogapsin1") == 0)
      {
        scoring.no_gaps_in_a = 1;
      }
      else if(strcasecmp(argv[argi], "--nogapsin2") == 0)
      {
        scoring.no_gaps_in_b = 1;
      }
      else if(strcasecmp(argv[argi], "--nomismatches") == 0)
      {
        scoring.no_mismatches = 1;
      }
      else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
      {
        // Already dealt with
        //case_sensitive = 1;
      }
      else if(strcasecmp(argv[argi], "--printscores") == 0)
      {
        print_scores = 1;
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
      else if(strcasecmp(argv[argi], "--zam") == 0)
      {
        print_zam = 1;
      }
      else if(strcasecmp(argv[argi], "--stdin") == 0)
      {
        // Similar to --file argument below
        cmdline_add_files("-", NULL);
      }
      else if(argi == argc-1)
      {
        // All the remaining options take an extra argument
        print_usage("Unknown argument without parameter: %s",argv[argi]);
      }
      else if(strcasecmp(argv[argi], "--scoring") == 0)
      {
        // This handled above
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_matrix") == 0)
      {
        gzFile sub_matrix_file = gzopen(argv[argi+1], "r");
        // gzbuffer(sub_matrix_file, 16384); // doesn't seem to work

        align_scoring_load_matrix(sub_matrix_file, argv[argi+1],
                                  &scoring, case_sensitive);

        //gzclose_r(sub_matrix_file); // doesn't seem to work
        gzclose(sub_matrix_file);
        substitutions_set = 1;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_pairs") == 0)
      {
        gzFile sub_pairs_file = gzopen(argv[argi+1], "r");
        //gzbuffer(sub_pairs_file, 16384); // doesn't seem to work
        
        align_scoring_load_pairwise(sub_pairs_file, argv[argi+1],
                                    &scoring, case_sensitive);
        
        //gzclose_r(sub_pairs_file); // doesn't seem to work
        gzclose(sub_pairs_file);
        substitutions_set = 1;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--match") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring.match))
        {
          print_usage("Invalid --match argument ('%s') must be an int",
                      argv[argi+1]);
        }

        match_set = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--mismatch") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring.mismatch))
        {
          print_usage("Invalid --mismatch argument ('%s') must be an int",
                      argv[argi+1]);
        }

        mismatch_set = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring.gap_open))
        {
          print_usage("Invalid --gapopen argument ('%s') must be an int",
                      argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring.gap_extend))
        {
          print_usage("Invalid --gapextend argument ('%s') must be an int",
                      argv[argi+1]);
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
        else if(strcmp(argv[argi+1], "-") == 0 && strcmp(argv[argi+2], "-") == 0)
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

        scoring_add_wildcard(&scoring, argv[argi+1][0], wildscore);

        argi += 2; // took two arguments
      }
      else
      {
        // Error - unknown option
        print_usage("Unknown argument '%s'", argv[argi]);
      }
    }
    else
    {
      if(argc - argi != 2)
      {
        print_usage("Unknown options: '%s'", argv[argi]);
      }
      break;
    }
  }

  if((match_set && !mismatch_set && !scoring.no_mismatches) ||
     (!match_set && mismatch_set))
  {
    print_usage("--match --mismatch must both be set or neither set");
  }
  else if(substitutions_set && !match_set)
  {
    // if substitution table set and not match/mismatch
    scoring.use_match_mismatch = 0;
  }

  if((scoring.no_gaps_in_a || scoring.no_gaps_in_b) && scoring.no_mismatches)
  {
    print_usage("--nogaps.. --nomismatches cannot be used at together");
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

  if(print_zam && (print_pretty || print_scores || print_colour || print_fasta))
  {
    print_usage("Cannot use --printscore, --printfasta, --pretty or --colour "
                "with --zam");
  }
  // End of set up

  // Align!
  nw = needleman_wunsch_new();
  result = alignment_create(256);

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

  // Free memory for storing alignment results
  needleman_wunsch_free(nw);
  alignment_free(result);

  return EXIT_SUCCESS;
}
