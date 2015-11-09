/*
 alignment_cmdline.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// Turn on debugging output by defining SEQ_ALIGN_VERBOSE
//#define SEQ_ALIGN_VERBOSE

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <limits.h> // INT_MIN
#include <stdarg.h> // for va_list

#include "seq_file/seq_file.h"

#include "alignment.h"
#include "alignment_cmdline.h"
#include "alignment_scoring_load.h"

// File loading
int file_list_length = 0;
int file_list_capacity = 0;
char **file_paths1 = NULL, **file_paths2 = NULL;

char parse_entire_int(char *str, int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (int)tmp;
    return 1;
  }
}

char parse_entire_uint(char *str, unsigned int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(tmp > UINT_MAX || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (unsigned int)tmp;
    return 1;
  }
}

static void print_usage(enum SeqAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...)
  __attribute__((format(printf, 4, 5)))
  __attribute__((noreturn));

static void print_usage(enum SeqAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fprintf(stderr, "\n");
  }

  fprintf(stderr, "usage: %s [OPTIONS] [seq1 seq2]\n", cmdstr);

  fprintf(stderr,
"  %s optimal %s alignment (maximises score).  \n"
"  Takes a pair of sequences on the command line, or can read from a\n"
"  file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.\n\n",
          cmd_type == SEQ_ALIGN_SW_CMD ? "Smith-Waterman" : "Needleman-Wunsch",
          cmd_type == SEQ_ALIGN_SW_CMD ? "local" : "global");

  fprintf(stderr,
"  OPTIONS:\n"
"    --file <file>        Sequence file reading with gzip support - read two\n"
"                         sequences at a time and align them\n"
"    --files <f1> <f2>    Read one sequence from each file to align at one time\n"
"    --stdin              Read from STDIN (same as '--file -')\n"
"\n"
"    --case_sensitive     Use case sensitive character comparison [default: off]\n"
"\n"
"    --match <score>      [default: %i]\n"
"    --mismatch <score>   [default: %i]\n"
"    --gapopen <score>    [default: %i]\n"
"    --gapextend <score>  [default: %i]\n"
"\n"
"    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>\n"
"    --substitution_matrix <file>  see details for formatting\n"
"    --substitution_pairs <file>   see details for formatting\n"
"\n"
"    --wildcard <w> <s>   Character <w> matches all characters with score <s>\n\n",
          defaults[0], defaults[1],
          defaults[2], defaults[3]);

  if(cmd_type == SEQ_ALIGN_SW_CMD)
  {
    // SW specific
    fprintf(stderr,
"    --minscore <score>   Minimum required score\n"
"                         [default: match * MAX(0.2 * length, 2)]\n"
"    --maxhits <hits>     Maximum number of results per alignment\n"
"                         [default: no limit]\n"
"\n"
"    --context <n>        Print <n> bases of context\n"
"    --printseq           Print sequences before local alignments\n");
  }
  else
  {
    // NW specific
    fprintf(stderr,
"\n"
"    --freestartgap       No penalty for gap at start of alignment\n"
"    --freeendgap         No penalty for gap at end of alignment\n"
"\n"
"    --printscores        Print optimal alignment scores\n"
"    --zam                A funky type of output\n");
  }

  fprintf(stderr,
"    --printmatrices      Print dynamic programming matrices\n"
"    --printfasta         Print fasta header lines\n"
"    --pretty             Print with a descriptor line\n"
"    --colour             Print with colour\n"
"\n"
"  Experimental Options:\n"
"    --nogapsin1          No gaps allowed within the first sequence\n"
"    --nogapsin2          No gaps allowed within the second sequence\n"
"    --nogaps             No gaps allowed in either sequence\n");

  fprintf(stderr,
"    --nomismatches       No mismatches allowed%s\n",
          cmd_type == SEQ_ALIGN_SW_CMD ? "" : " (cannot be used with --nogaps..)");

  printf(
"\n"
" DETAILS:\n"
"  * For help choosing scoring, see the README file. \n"
"  * Gap (of length N) penalty is: (open+N*extend)\n"
"  * To do alignment without affine gap penalty, set '--gapopen 0'.\n"
"  * Scoring files should be matrices, with entries separated by a single\n"
"    character or whitespace. See files in the 'scores' directory for examples.\n"
"\n"
"  turner.isaac@gmail.com  (compiled: %s %s)\n", __DATE__, __TIME__);

  exit(EXIT_FAILURE);
}

void cmdline_free(cmdline_t *cmd)
{
  free(cmd->file_paths1);
  free(cmd->file_paths2);
  free(cmd);
}

#define usage(fmt,...) print_usage(cmd_type,defaults,argv[0],fmt, ##__VA_ARGS__)

cmdline_t* cmdline_new(int argc, char **argv, scoring_t *scoring,
                       enum SeqAlignCmdType cmd_type)
{
  cmdline_t* cmd = calloc(1, sizeof(cmdline_t));
  cmd->file_list_length = 0;
  cmd->file_list_capacity = 256;
  cmd->file_paths1 = malloc(sizeof(char*) * cmd->file_list_capacity);
  cmd->file_paths2 = malloc(sizeof(char*) * cmd->file_list_capacity);
  cmd->seq1 = cmd->seq2 = NULL;
  // All values initially 0

  // Store defaults
  score_t defaults[4] = {scoring->match, scoring->mismatch,
                         scoring->gap_open, scoring->gap_extend};

  if(argc == 1) usage(NULL);

  // First run through arguments to set up case_sensitive and scoring system

  // case sensitive needs to be dealt with first
  // (is is used to construct hash table for swap_scores)
  char scoring_set = 0, substitutions_set = 0, match_set = 0, mismatch_set = 0;

  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(strcasecmp(argv[argi], "--help") == 0 ||
       strcasecmp(argv[argi], "-help") == 0 ||
       strcasecmp(argv[argi], "-h") == 0)
    {
      usage(NULL);
    }
    else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
    {
      cmd->case_sensitive = 1;
    }
    else if(strcasecmp(argv[argi], "--scoring") == 0)
    {
      if(scoring_set)
      {
        usage("More than one scoring system specified - not permitted");
      }

      if(strcasecmp(argv[argi+1], "PAM30") == 0)
      {
        scoring_system_PAM30(scoring);
      }
      else if(strcasecmp(argv[argi+1], "PAM70") == 0)
      {
        scoring_system_PAM70(scoring);
      }
      else if(strcasecmp(argv[argi+1], "BLOSUM80") == 0)
      {
        scoring_system_BLOSUM80(scoring);
      }
      else if(strcasecmp(argv[argi+1], "BLOSUM62") == 0)
      {
        scoring_system_BLOSUM62(scoring);
      }
      else if(strcasecmp(argv[argi+1], "DNA_HYBRIDIZATION") == 0)
      {
        scoring_system_DNA_hybridization(scoring);
      }
      else {
        usage("Unknown --scoring choice, not one of "
              "PAM30|PAM70|BLOSUM80|BLOSUM62");
      }

      scoring_set = 1;
      argi++; // took an argument
    }
  }

  for(argi = 1; argi < argc; argi++)
  {
    if(argv[argi][0] == '-')
    {
      // strcasecmp does case insensitive comparison
      if(strcasecmp(argv[argi], "--freestartgap") == 0)
      {
        if(cmd_type != SEQ_ALIGN_NW_CMD)
          usage("--freestartgap only valid with Needleman-Wunsch");
        scoring->no_start_gap_penalty = true;
      }
      else if(strcasecmp(argv[argi], "--freeendgap") == 0)
      {
        if(cmd_type != SEQ_ALIGN_NW_CMD)
          usage("--freeendgap only valid with Needleman-Wunsch");
        scoring->no_end_gap_penalty = true;
      }
      else if(strcasecmp(argv[argi], "--nogaps") == 0)
      {
        scoring->no_gaps_in_a = true;
        scoring->no_gaps_in_b = true;
      }
      else if(strcasecmp(argv[argi], "--nogapsin1") == 0)
      {
        scoring->no_gaps_in_a = true;
      }
      else if(strcasecmp(argv[argi], "--nogapsin2") == 0)
      {
        scoring->no_gaps_in_b = true;
      }
      else if(strcasecmp(argv[argi], "--nomismatches") == 0)
      {
        scoring->no_mismatches = true;
      }
      else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
      {
        // Already dealt with
        //case_sensitive = true;
      }
      else if(strcasecmp(argv[argi], "--printseq") == 0)
      {
        if(cmd_type != SEQ_ALIGN_SW_CMD)
          usage("--printseq only valid with Smith-Waterman");
        cmd->print_seq = true;
      }
      else if(strcasecmp(argv[argi], "--printmatrices") == 0)
      {
        cmd->print_matrices = true;
      }
      else if(strcasecmp(argv[argi], "--printscores") == 0)
      {
        if(cmd_type != SEQ_ALIGN_NW_CMD)
          usage("--printscores only valid with Needleman-Wunsch");
        cmd->print_scores = true;
      }
      else if(strcasecmp(argv[argi], "--printfasta") == 0)
      {
        cmd->print_fasta = true;
      }
      else if(strcasecmp(argv[argi], "--pretty") == 0)
      {
        cmd->print_pretty = true;
      }
      else if(strcasecmp(argv[argi], "--colour") == 0)
      {
        cmd->print_colour = true;
      }
      else if(strcasecmp(argv[argi], "--zam") == 0)
      {
        if(cmd_type != SEQ_ALIGN_NW_CMD)
          usage("--zam only valid with Needleman-Wunsch");
        cmd->zam_stle_output = true;
      }
      else if(strcasecmp(argv[argi], "--stdin") == 0)
      {
        // Similar to --file argument below
        cmdline_add_files(cmd, "", NULL);
        cmd->interactive = true;
      }
      else if(argi == argc-1)
      {
        // All the remaining options take an extra argument
        usage("Unknown argument without parameter: %s", argv[argi]);
      }
      else if(strcasecmp(argv[argi], "--scoring") == 0)
      {
        // This handled above
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_matrix") == 0)
      {
        gzFile sub_matrix_file = gzopen(argv[argi+1], "r");
        if(sub_matrix_file == NULL) usage("Couldn't read: %s", argv[argi+1]);

        align_scoring_load_matrix(sub_matrix_file, argv[argi+1],
                                  scoring, cmd->case_sensitive);

        gzclose(sub_matrix_file);
        substitutions_set = true;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_pairs") == 0)
      {
        gzFile sub_pairs_file = gzopen(argv[argi+1], "r");
        if(sub_pairs_file == NULL) usage("Couldn't read: %s", argv[argi+1]);

        align_scoring_load_pairwise(sub_pairs_file, argv[argi+1],
                                    scoring, cmd->case_sensitive);

        gzclose(sub_pairs_file);
        substitutions_set = true;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--minscore") == 0)
      {
        if(cmd_type != SEQ_ALIGN_SW_CMD)
          usage("--minscore only valid with Smith-Waterman");

        if(!parse_entire_int(argv[argi+1], &cmd->min_score))
          usage("Invalid --minscore <score> argument (must be a +ve int)");

        cmd->min_score_set = true;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--maxhits") == 0)
      {
        if(cmd_type != SEQ_ALIGN_SW_CMD)
          usage("--maxhits only valid with Smith-Waterman");

        if(!parse_entire_uint(argv[argi+1], &cmd->max_hits_per_alignment))
          usage("Invalid --maxhits <hits> argument (must be a +ve int)");

        cmd->max_hits_per_alignment_set = true;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--context") == 0)
      {
        if(cmd_type != SEQ_ALIGN_SW_CMD)
          usage("--context only valid with Smith-Waterman");

        if(!parse_entire_uint(argv[argi+1], &cmd->print_context))
          usage("Invalid --context <c> argument (must be >= 0)");

        argi++;
      }
      else if(strcasecmp(argv[argi], "--match") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->match))
        {
          usage("Invalid --match argument ('%s') must be an int", argv[argi+1]);
        }

        match_set = true;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--mismatch") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->mismatch))
        {
          usage("Invalid --mismatch argument ('%s') must be an int", argv[argi+1]);
        }

        mismatch_set = true;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_open))
        {
          usage("Invalid --gapopen argument ('%s') must be an int", argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_extend))
        {
          usage("Invalid --gapextend argument ('%s') must be an int",
                argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--file") == 0)
      {
        cmdline_add_files(cmd, argv[argi+1], NULL);
        argi++; // took an argument
      }
      // Remaining options take two arguments but check themselves
      else if(strcasecmp(argv[argi], "--files") == 0)
      {
        if(argi >= argc-2)
        {
          usage("--files option takes 2 arguments");
        }
        else if(strcmp(argv[argi+1], "-") == 0 && strcmp(argv[argi+2], "-") == 0)
        {
          // Read both from stdin
          cmdline_add_files(cmd, argv[argi+1], NULL);
        }
        else
        {
          cmdline_add_files(cmd, argv[argi+1], argv[argi+2]);
        }

        argi += 2; // took two arguments
      }
      else if(strcasecmp(argv[argi], "--wildcard") == 0)
      {
        int wildscore = 0;

        if(argi == argc-2 || strlen(argv[argi+1]) != 1 ||
           !parse_entire_int(argv[argi+2], &wildscore))
        {
          usage("--wildcard <w> <s> takes a single character and a number");
        }

        scoring_add_wildcard(scoring, argv[argi+1][0], wildscore);

        argi += 2; // took two arguments
      }
      else usage("Unknown argument '%s'", argv[argi]);
    }
    else
    {
      if(argc - argi != 2) usage("Unknown options: '%s'", argv[argi]);
      break;
    }
  }

  if((match_set && !mismatch_set && !scoring->no_mismatches) ||
     (!match_set && mismatch_set))
  {
    usage("--match --mismatch must both be set or neither set");
  }
  else if(substitutions_set && !match_set)
  {
    // if substitution table set and not match/mismatch
    scoring->use_match_mismatch = 0;
  }

  if(scoring->use_match_mismatch && scoring->match < scoring->mismatch) {
    usage("Match value should not be less than mismatch penalty");
  }

  // Cannot guarantee that we can perform a global alignment if nomismatches
  // and nogaps is true
  if(cmd_type == SEQ_ALIGN_NW_CMD && scoring->no_mismatches &&
     (scoring->no_gaps_in_a || scoring->no_gaps_in_b))
  {
    usage("--nogaps.. --nomismatches cannot be used at together");
  }

  // Check for extra unused arguments
  // and set seq1 and seq2 if they have been passed
  if(argi < argc)
  {
    cmd->seq1 = argv[argi];
    cmd->seq2 = argv[argi+1];
  }

  if(cmd->seq1 == NULL && cmd->file_list_length == 0)
  {
    usage("No input specified");
  }

  if(cmd->zam_stle_output &&
     (cmd->print_pretty || cmd->print_scores ||
      cmd->print_colour || cmd->print_fasta))
  {
    usage("Cannot use --printscore, --printfasta, --pretty or --colour with "
          "--zam");
  }

  return cmd;
}


void cmdline_add_files(cmdline_t *cmd, char* p1, char* p2)
{
  if(cmd->file_list_length == cmd->file_list_capacity)
  {
    cmd->file_list_capacity *= 2;
    size_t mem = sizeof(char*) * cmd->file_list_capacity;
    cmd->file_paths1 = realloc(cmd->file_paths1, mem);
    cmd->file_paths2 = realloc(cmd->file_paths2, mem);

    if(cmd->file_paths1 == NULL || cmd->file_paths2 == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  cmd->file_paths1[cmd->file_list_length] = p1;
  cmd->file_paths2[cmd->file_list_length] = p2;
  cmd->file_list_length++;
}

size_t cmdline_get_num_of_file_pairs(cmdline_t *cmd)
{
  return cmd->file_list_length;
}

char* cmdline_get_file1(cmdline_t *cmd, size_t i)
{
  return cmd->file_paths1[i];
}

char* cmdline_get_file2(cmdline_t *cmd, size_t i)
{
  return cmd->file_paths2[i];
}

static seq_file_t* open_seq_file(const char *path, bool use_zlib)
{
  return (strcmp(path,"-") != 0 || use_zlib) ? seq_open(path)
                                             : seq_dopen(fileno(stdin), false, false, 0);
}

// If seq2 is NULL, read pair of entries from first file
// Otherwise read an entry from each
void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2),
                     bool use_zlib)
{
  seq_file_t *sf1, *sf2;

  if((sf1 = open_seq_file(path1, use_zlib)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path1);
    fflush(stderr);
    return;
  }

  if(path2 == NULL)
  {
    sf2 = sf1;
  }
  else if((sf2 = open_seq_file(path2, use_zlib)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path1);
    fflush(stderr);
    return;
  }

  // fprintf(stderr, "File buffer %zu zlib: %i\n", sf1->in.size, seq_use_gzip(sf1));

  read_t read1, read2;
  seq_read_alloc(&read1);
  seq_read_alloc(&read2);

  // Loop while we can read a sequence from the first file
  unsigned long alignments;

  for(alignments = 0; seq_read(sf1, &read1) > 0; alignments++)
  {
    if(seq_read(sf2, &read2) <= 0)
    {
      fprintf(stderr, "Alignment Error: Odd number of sequences - "
                      "I read in pairs!\n");
      fflush(stderr);
      break;
    }

    (align)(&read1, &read2);
  }

  // warn if no bases read
  if(alignments == 0)
  {
    fprintf(stderr, "Alignment Warning: empty input\n");
    fflush(stderr);
  }

  // Close files
  seq_close(sf1);

  if(path2 != NULL)
    seq_close(sf2);

  // Free memory
  seq_read_dealloc(&read1);
  seq_read_dealloc(&read2);
}
