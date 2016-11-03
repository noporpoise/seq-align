seq-align
=========
Smith-Waterman & Needleman-Wunsch Alignment in C  
url: https://github.com/noporpoise/seq-align  
author: Isaac Turner <turner.isaac@gmail.com>  
license: Public Domain  
updated: 18 Aug 2015

[![Build Status](https://travis-ci.org/noporpoise/seq-align.png?branch=master)](https://travis-ci.org/noporpoise/seq-align)

About
=====

C implementations of optimal local (Smith-Waterman) and global
(Needleman-Wunsch) alignment algorithms.  Written to be fast, portable and easy
to use.  Commandline utilities `smith_waterman` and `needleman_wunsch` provide
great flexibility.  Code can also be included easily in third party programs,
see `nw_example/` and `sw_example/` directories for examples.  Also includes perl
modules to provide perl API to programs.  

Features:
* Align any pair of ASCII sequences (DNA, Protein, words, etc.)
* Specify alignment scoring systems or choose a common one (BLOSUM etc.)
* Specify wildcards that match any character with a given score
* Align with no mismatches (`--nomismatches`), or no gaps (`--nogaps`) with
  local and global alignment.  When both used for local alignment, lists common
  substrings in order of length/score
* Read fastq, fasta, sam, bam, plain (one sequence per line) and gzipped files
* Display in colour (`--colour`)
* Show alignment context when doing local alignment (`--context <n>`)
* Allow a penalty free gap at the beginning or end of a global alignment
  (`--freestartgap` and `--freeendgap`)
* Default behaviour is case-insensitve matching, change using `--case_sensitive`

Build
-----

Download and build seq-align:

    $ git clone --recursive https://github.com/noporpoise/seq-align
    $ make

To update:

    $ git pull
    $ git submodule update --init

To run tests:

    $ make test

For those interested, the bundled depedencies used are:

* string_buffer [https://github.com/noporpoise/string_buffer]
* seq_file [https://github.com/noporpoise/seq_file]
* sort_r [https://github.com/noporpoise/sort_r]

To include seq-align in your own applications, look at the examples in `./examples/`. Perl wrapper modules are also available in `./perl/`.

Examples
========

Baiscs:

    $ ./needleman_wunsch CAGACGT CGATA
    C-AGACGT
    CGATA---

Print alignment scores:

    $ ./needleman_wunsch --printscores CAGACGT CGATA
    C-AGACGT
    CGATA---
    score: -11

Read from file (dna.fa.gz):

    >seqA
    ACAATAGAC
    >seqB
    ACGAATAGAT
    >seqC
    ACGTGA
    CAGAT
    >seqD
    GTGGACG
    AGTA

    $ ./needleman_wunsch --printscores --file dna.fa.gz
    AC-AATAGAC
    ACGAATAGAT
    score: -3

    ACGTGACAGAT
    GTGGACGAGTA
    score: -5


Reading from STDIN:

    $ gzip -d -c dna.fa.gz | ./needleman_wunsch --stdin
    AC-AATAGAC
    ACGAATAGAT

    ACGTGACAGAT
    GTGGACGAGTA

is the same as:

    $ gzip -d -c dna.fa.gz | ./needleman_wunsch --file -

Set different scoring systems:

    $ ./needleman_wunsch --match 1 --mismatch 0 --gapopen -10 --gapextend 0 ACGTGCCCCACAGAT AGGTGGACGAGAT


Smith-Waterman
==============

        usage: ./bin/smith_waterman [OPTIONS] [seq1 seq2]
          Smith-Waterman optimal local alignment (maximises score).  
          Takes a pair of sequences on the command line, or can read from a
          file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.

          OPTIONS:
            --file <file>        Sequence file reading with gzip support - read two
                                 sequences at a time and align them
            --files <f1> <f2>    Read one sequence from each file to align at one time
            --stdin              Read from STDIN (same as '--file -')

            --case_sensitive     Use case sensitive character comparison [default: off]

            --match <score>      [default: 2]
            --mismatch <score>   [default: -2]
            --gapopen <score>    [default: -1]
            --gapextend <score>  [default: -1]

            --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>
            --substitution_matrix <file>  see details for formatting
            --substitution_pairs <file>   see details for formatting

            --wildcard <w> <s>   Character <w> matches all characters with score <s>

            --minscore <score>   Minimum required score
                                 [default: match * MAX(0.2 * length, 2)]
            --maxhits <hits>     Maximum number of results per alignment
                                 [default: no limit]

            --context <n>        Print <n> bases of context
            --printseq           Print sequences before local alignments
            --printmatrices      Print dynamic programming matrices
            --printfasta         Print fasta header lines
            --pretty             Print with a descriptor line
            --colour             Print with colour

          Experimental Options:
            --nogapsin1          No gaps allowed within the first sequence
            --nogapsin2          No gaps allowed within the second sequence
            --nogaps             No gaps allowed in either sequence
            --nomismatches       No mismatches allowed

         DETAILS:
          * For help choosing scoring, see the README file. 
          * Gap (of length N) penalty is: (open+N*extend)
          * To do alignment without affine gap penalty, set '--gapopen 0'.
          * Scoring files should be matrices, with entries separated by a single
            character or whitespace. See files in the 'scores' directory for examples.

          turner.isaac@gmail.com  (compiled: Feb 11 2015 20:20:34)


Needleman-Wunsch
================

        usage: ./bin/needleman_wunsch [OPTIONS] [seq1 seq2]
          Needleman-Wunsch optimal global alignment (maximises score).  
          Takes a pair of sequences on the command line, or can read from a
          file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.

          OPTIONS:
            --file <file>        Sequence file reading with gzip support - read two
                                 sequences at a time and align them
            --files <f1> <f2>    Read one sequence from each file to align at one time
            --stdin              Read from STDIN (same as '--file -')

            --case_sensitive     Use case sensitive character comparison [default: off]

            --match <score>      [default: 1]
            --mismatch <score>   [default: -2]
            --gapopen <score>    [default: -4]
            --gapextend <score>  [default: -1]

            --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>
            --substitution_matrix <file>  see details for formatting
            --substitution_pairs <file>   see details for formatting

            --wildcard <w> <s>   Character <w> matches all characters with score <s>


            --freestartgap       No penalty for gap at start of alignment
            --freeendgap         No penalty for gap at end of alignment

            --printscores        Print optimal alignment scores
            --zam                A funky type of output
            --printmatrices      Print dynamic programming matrices
            --printfasta         Print fasta header lines
            --pretty             Print with a descriptor line
            --colour             Print with colour

          Experimental Options:
            --nogapsin1          No gaps allowed within the first sequence
            --nogapsin2          No gaps allowed within the second sequence
            --nogaps             No gaps allowed in either sequence
            --nomismatches       No mismatches allowed (cannot be used with --nogaps..)

         DETAILS:
          * For help choosing scoring, see the README file. 
          * Gap (of length N) penalty is: (open+N*extend)
          * To do alignment without affine gap penalty, set '--gapopen 0'.
          * Scoring files should be matrices, with entries separated by a single
            character or whitespace. See files in the 'scores' directory for examples.

          turner.isaac@gmail.com  (compiled: Feb 11 2015 20:20:34)


Longest Common Substring
========================

    ./bin/lcs [options] <sequence>
      Print substrings in order of length

Print maximal substrings from longest to shortest


Scoring Penalties
-----------------

Proteins:

| Query Length  | Substitution Matrix | Gap Costs (open,extend) |
|---------------|---------------------|-------------------------|
| <35           | PAM-30              | (9,1)                   |
| 35-50         | PAM-70              | (10,1)                  |
| 50-85         | BLOSUM-80           | (10,1)                  |
| 85            | BLOSUM-62           | (10,1)                  |

[Table from: http://www.ncbi.nlm.nih.gov/blast/html/sub_matrix.html]

    gap (of length N) penalty: gap_open + N*gap_extend

NCBI BLAST Quote [from: http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml#Reward-penalty]:

    Many nucleotide searches use a simple scoring system that consists of a "reward"
    for a match and a "penalty" for a mismatch. The (absolute) reward/penalty ratio
    should be increased as one looks at more divergent sequences. A ratio of 0.33
    (1/-3) is appropriate for sequences that are about 99% conserved; a ratio of 0.5
    (1/-2) is best for sequences that are 95% conserved; a ratio of about one (1/-1)
    is best for sequences that are 75% conserved [1].

NCBI Gap (open, extend) values:

| open | extend |
|------|--------|
| -5   | -2     |
| -2   | -2     |
| -1   | -2     |
| -0   | -2     |
| -3   | -1     |
| -2   | -1     |
| -1   | -1     |

Our default (for now) are:  
gap_open/gap_extend: (-4,-1)  
match/mismatch: (1,-2)

 for help on choosing scoring matrices see:
 http://www.ebi.ac.uk/help/matrix.html

License
=======

 NCBI protein align matices from:
 ftp://ftp.ncbi.nih.gov/blast/matrices/

seq-align: Public Domain
    This is free and unencumbered software released into the public domain.

    Anyone is free to copy, modify, publish, use, compile, sell, or
    distribute this software, either in source code form or as a compiled
    binary, for any purpose, commercial or non-commercial, and by any
    means.

    In jurisdictions that recognize copyright laws, the author or authors
    of this software dedicate any and all copyright interest in the
    software to the public domain. We make this dedication for the benefit
    of the public at large and to the detriment of our heirs and
    successors. We intend this dedication to be an overt act of
    relinquishment in perpetuity of all present and future rights to this
    software under copyright law.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
    OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

    For more information, please refer to <http://unlicense.org>


DEVELOPMENT
===========

Feel free to contact me to request features.  Bug reports are appreciated.  
(turner.isaac@gmail.com)

Possible features:
1. May add option to have no gaps in alignment result for longest sequence
   with `--nogapsinlongest`

Algorithms to investigate:
Gotoh
Hirschberg's algorithm using linear space
Myer's bit-vector algorithm
Combination of Myer's and Hirschberg's algorithm
