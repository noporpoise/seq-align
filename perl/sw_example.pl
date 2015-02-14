#!/usr/bin/env perl

# perl/sw_example.pl
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

use SmithWaterman;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./sw_example.pl <seq1> <seq2>\n";

  exit;
}

my $sw = new SmithWaterman('minscore' => 2);

my @aligns = (['asdfasdf','asdfasdf'],
              ['dogg', 'ggod'],
              ['asdfas', 'fadsas'],
              ['hello henry', 'wellohenry'],
              ['asdf', 'lkj']);

for my $align (@aligns)
{
  print "Doing '$align->[0]' vs '$align->[1]'\n";
  $sw->do_alignment(@$align);

  #while
  if(defined(my $hit = $sw->get_next_hit()))
  {
    print "hit:\n";
    print "  ".$hit->{'align1'}."\n";
    print "  ".$hit->{'sep'}."\n";
    print "  ".$hit->{'align2'}."\n";
  }
}

$sw->destructor();
