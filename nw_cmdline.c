/*
 nw_cmdline.c
 project: NeedlemanWunsch
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/needlemanwunsch
 Copyright (C) 06-Dec-2011
 
 see: README

 == License
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower
#include <stdarg.h> // required for va_list

// my utility functions
#include "string_buffer.h"
#include "seq_file.h"
#include "utility_lib.h"

// Alignment scoring and loading
#include "alignment_scoring_load.h"
#include "alignment_cmdline.h"

#include "needleman_wunsch.h"

// For this run
char* cmd;
char print_colour = 0, print_pretty = 0, print_scores = 0,
     print_fasta = 0, print_zam = 0;

SCORING_SYSTEM* scoring = NULL;

// Alignment results stored here
char *alignment_a = NULL, *alignment_b = NULL;
t_buf_pos alignment_max_length;

void set_default_scoring()
{
  scoring = scoring_system_default();
}

void print_usage(char* err_fmt, ...)
{
  if(err_fmt != NULL)
  {
    va_list argptr;
    va_start(argptr, err_fmt);
    
    StrBuf *error = strbuf_init(200);
    strbuf_append_str(error, "NeedlemanWunsch Error: ");
    strbuf_vsprintf(error, strbuf_len(error), err_fmt, argptr);
    strbuf_chomp(error);

    va_end(argptr);

    fprintf(stderr, "%s\n", error->buff);
    fprintf(stderr, "Use -h option to print help\n");
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
"    --nogaps             No gaps allowed in the alignment\n"
"    --nomismatches       No mismatches allowed - not to be used with --nogaps\n"
"\n"
"    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>\n"
"    --substitution_matrix <file>  see details for formatting\n"
"    --substitution_pairs <file>   see details for formatting\n"
"\n"
"    --wildcard <w> <s>   Character <w> matches all characters with score <s>\n",
          scoring->match, scoring->mismatch,
          scoring->gap_open, scoring->gap_extend);

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

void align_zam(const char *seq_a, const char *seq_b)
{
  needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

  // Swap '-' for '_'
  int i;
  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-')
    {
      alignment_a[i] = '_';
    }

    if(alignment_b[i] == '-')
    {
      alignment_b[i] = '_';
    }
  }

  int num_of_mismatches = 0;
  int num_of_indels = 0;

  // Print branch 1
  printf("Br1:%s\n", alignment_a);

  // Print spacer
  printf("    ");

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '_' || alignment_b[i] == '_')
    {
      printf(" ");
      num_of_indels++;
    }
    else if((scoring->case_sensitive && alignment_a[i] != alignment_b[i]) ||
            tolower(alignment_a[i]) != tolower(alignment_b[i]))
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
  printf("Br2:%s\n", alignment_b);

  // print mismatch indel numbers
  printf("%i %i\n\n", num_of_mismatches, num_of_indels);
}

void align(const char *seq_a, const char *seq_b,
           const char *seq_a_name, const char *seq_b_name)
{
  if(print_zam)
  {
    align_zam(seq_a, seq_b);
    fflush(stdout);
    return;
  }

  int score = needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

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
    alignment_colour_print_against(alignment_a, alignment_b,
                                   scoring->case_sensitive);
  }
  else
  {
    printf("%s", alignment_a);
  }
  printf("\n");
  
  if(print_pretty)
  {
    // Print spacer
    alignment_print_spacer(alignment_a, alignment_b, scoring);

    printf("\n");
  }
  else if(print_fasta && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_colour)
  {
    // Print alignment line 2
    alignment_colour_print_against(alignment_b, alignment_a,
                                   scoring->case_sensitive);
  }
  else
  {
    printf("%s", alignment_b);
  }
  printf("\n");

  if(print_scores)
  {
    printf("score: %i\n", score);
  }
  
  printf("\n");

  fflush(stdout);
}

void align_pair_from_file(StrBuf *seq1, StrBuf *seq2,
                          const char *seq1_name, const char *seq2_name)
{
  // Check memory
  t_buf_pos new_max_alignment = seq1->len + seq2->len;

  if(new_max_alignment > alignment_max_length)
  {
    // Expand memory used for storing result
    alignment_max_length = new_max_alignment;

    if(!nw_realloc_mem((unsigned int)new_max_alignment,
                       &alignment_a, &alignment_b))
    {
      print_usage("Ran out of memory");
    }
  }

  align(seq1->buff, seq2->buff, seq1_name, seq2_name);
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
  // (is is used to construct hash table for swap_table)
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
      if(strcasecmp(argv[argi], "--freestartgap") == 0)
      {
        scoring->no_start_gap_penalty = 1;
      }
      else if(strcasecmp(argv[argi], "--freeendgap") == 0)
      {
        scoring->no_end_gap_penalty = 1;
      }
      else if(strcasecmp(argv[argi], "--nogaps") == 0)
      {
        scoring->no_gaps = 1;
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
          print_usage("Invalid --mismatch argument ('%s') must be an int",
                      argv[argi+1]);
        }

        mismatch_set = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_open))
        {
          print_usage("Invalid --gapopen argument ('%s') must be an int",
                      argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_extend))
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

        scoring_add_wildcard(scoring, argv[argi+1][0], wildscore);

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

  if(match_set != mismatch_set)
  {
    print_usage("--match --mismatch must both be set or neither set");
  }
  else if(substitutions_set && !match_set)
  {
    // if substitution table set and not match/mismatch
    scoring->use_match_mismatch = 0;
  }

  if(scoring->no_gaps && scoring->no_mismatches)
  {
    print_usage("--nogaps --nomismatches cannot be used at together");
  }

  if(scoring->no_gaps)
  {
    scoring->no_start_gap_penalty = 1;
    scoring->no_end_gap_penalty = 1;
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
  if(seq1 != NULL)
  {
    // Align seq1 and seq2
    alignment_max_length = nw_alloc_mem(seq1, seq2, &alignment_a, &alignment_b);
    align(seq1, seq2, NULL, NULL);
  }
  else
  {
    // Set up default memory for aligning from stdin / files
    alignment_max_length = 1000; // equivalent to two strings of 500bp
    alignment_a = (char*) malloc((alignment_max_length+1) * sizeof(char));
    alignment_b = (char*) malloc((alignment_max_length+1) * sizeof(char));
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
  free(alignment_a);
  free(alignment_b);

  scoring_free(scoring);

  return EXIT_SUCCESS;
}
