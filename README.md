seq-align
=========
Smith-Waterman & Needleman-Wunsch Alignment in C    
url: https://github.com/noporpoise/seq-align  
author: Isaac Turner <turner.isaac@gmail.com>  
license: GPLv3
updated: 8 Oct 2012

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

Fetch and build dependencies (requires `git`, `make` and `zlib`)

    $ cd libs
    $ make
    $ cd ..

Build seq-align:

    $ make

For those interested, the depedencies used are:

* htslib [https://github.com/samtools/htslib]
* bit_array [https://github.com/noporpoise/BitArray]
* string_buffer [https://github.com/noporpoise/string_buffer]
* seq_file [https://github.com/noporpoise/seq_file]

Smith-Waterman
==============

    usage: ./smith_waterman [OPTIONS] [seq1 seq2]
      Smith-Waterman optimal local alignment (maximises score).  
      Takes a pair of sequences on the command line, or can read from a
      file and from sequence piped in.  Can read gzip files and FASTA.
    
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
    
        --nogaps             No gaps allowed in the alignment
        --nomismatches       No mismatches allowed. 
          When used together, prints longest common substrings in order of length
    
        --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>
        --substitution_matrix <file>  see details for formatting
        --substitution_pairs <file>   see details for formatting
    
        --wildcard <w> <s>   Character <w> matches all characters with score <s>
    
        --minscore <score>   Minimum required score
                             [default: match * MAX(0.2 * length, 2)]
        --maxhits <hits>     Maximum number of results per alignment
                             [default: no limit]
        --context <n>        Print <n> bases of context
        --printfasta         Print fasta header lines
        --printseq           Print sequences before local alignments
        --pretty             Print with a descriptor line
        --colour             Print with colour

Needleman-Wunsch
================

    usage: ./needleman_wunsch [OPTIONS] [seq1 seq2]
      Needleman-Wunsch optimal global alignment (maximises score). Takes a pair 
      of sequences on the command line, reads from a file and from sequence 
      piped in.  Can read gzip files and those in FASTA, FASTQ or plain format.
    
      OPTIONS:
        --file <file>        Sequence file reading with gzip support
        --files <f1> <f2>    Read one sequence from each file at a time to align
        --stdin              Read from STDIN (same as '--file -')
    
        --case_sensitive     Case sensitive character comparison
    
        --match <score>      [default: 1]
        --mismatch <score>   [default: -2]
        --gapopen <score>    [default: -4]
        --gapextend <score>  [default: -1]
    
        --nogaps             No gaps allowed in the alignment
        --nomismatches       No mismatches allowed - not to be used with --nogaps
    
        --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>
        --substitution_matrix <file>  see details for formatting
        --substitution_pairs <file>   see details for formatting
    
        --wildcard <w> <s>   Character <w> matches all characters with score <s>
    
        --freestartgap       No penalty for gap at start of alignment
        --freeendgap         No penalty for gap at end of alignment
    
        --printscores        Print optimal alignment scores
        --printfasta         Print fasta header lines
        --pretty             Print with a descriptor line
        --colour             Print with colour
        --zam                A funky type of output


Baiscs:

    $ ./needleman_wunsch CAGACGT CGATA
    CAGACGT
    C--GATA

Print alignment scores:

    $ ./needleman_wunsch --printscores CAGACGT CGATA
    CAGACGT
    C--GATA
    score: -15

Read from file (dna.fa):

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
-------

 NCBI protein align matices from:
 ftp://ftp.ncbi.nih.gov/blast/matrices/

GPLv3

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

DEVELOPMENT
===========

Feel free to contact me to request features.  Bug reports are appreciated.  
(turner.isaac@gmail.com)

1) Add tests
2) Finish adding the following options:

* No unneeded gaps:
--nogapsin1  
--nogapsin2  
--nogaps

* Penalised gaps at ends:
--gapsonlyatendsin1  
--gapsonlyatendsin2  
--gapsonlyatends

* free gaps at ends:
--freestartgap  
--freeendgap
