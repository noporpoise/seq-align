#!/usr/bin/env perl

# perl/nw_example.pl
# url: https://github.com/noporpoise/seq-align
# maintainer: Isaac Turner <turner.isaac@gmail.com>
# license: Public Domain, no warranty
# date: Oct 2012

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use NeedlemanWunsch;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./nw_example.pl <seq1> <seq2>\n";

  exit;
}

my $nw = new NeedlemanWunsch();

my @aligns = (['asdfasdf','asdfasdf'],
              ['dogg', 'ggod'],
              ['asdfas', 'fadsas'],
              ['', ''],
              ['asdf', 'lkj']);

for my $align (@aligns)
{
  print "Doing '$align->[0]' vs '$align->[1]'\n";
  my $result = $nw->do_alignment(@$align);
  $nw->print_alignment($result);
}

$nw->destructor();
