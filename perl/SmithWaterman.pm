package SmithWaterman;

# perl/SmithWaterman.pm
# url: https://github.com/noporpoise/seq-align
# maintainer: Isaac Turner <turner.isaac@gmail.com>
# license: Public Domain, no warranty
# date: Oct 2012

use strict;
use warnings;

use File::Basename; # dirname
use Carp; # for reporting warnings and errors
use FileHandle; # provides autoflush
use IPC::Open2; # for opening handle to a process

use constant {PROMPT_LINE => 'next [h]it or [a]lignment: '};

sub new
{
  my ($class, @args) = @_;

  my %options = (@args);

  my $cmdline = dirname(__FILE__)."/../bin/smith_waterman";
  my $timeout = 5;

  for my $key (keys %options)
  {
    $options{lc($key)} = $options{$key};
  }

  if(!defined($options{'nomismatches'}) &&
     defined($options{'match'}) != defined($options{'mismatch'}))
  {
    carp("Cannot set only one of match/mismatch");
  }

  if(defined($options{'cmd'}))
  {
    $cmdline = $options{'cmd'};
  }

  $cmdline .= " --stdin --pretty";

  while(my ($key,$value) = each(%options))
  {
    if(grep(/^$key$/i, qw(case_sensitive nogaps nomismatches)))
    {
      # 'Flag' options -- they have no args
      if($value)
      {
        $cmdline .= " --$key";
      }
    }
    elsif($key =~ /^scoring$/i)
    {
      if(!grep($value, qw(PAM30 PAM70 BLOSUM80 BLOSUM62)))
      {
        carp("Invalid SmithWaterman.pm scoring '$value'");
      }
      else
      {
        $cmdline .= " --scoring $value";
      }
    }
    elsif($key eq "timeout") {
      $timeout = $value;
    }
    elsif(grep(/^$key$/i, qw(substitution_matrix substitution_pairs
                             match mismatch gapopen gapextend
                             minscore maxhits wildcard)))
    {
      $cmdline .= " --$key $value";
    }
    else
    {
      carp("Unknown option '$key' => '$value' (ignored)");
    }
  }

  $cmdline .= " 2>&1";

  my ($in, $out, $err);

  #print "running '$cmdline'\n";

  my $pid = open2($in, $out, $cmdline)
    or die("Cannot run cmd: '$cmdline'");

  $in->autoflush();
  $out->autoflush();

  my $self = {
    _in => $in,
    _out => $out,
    _err => $err,
    _pid => $pid,
    _timeout => $timeout,
    _align_number => -1,
    _seq1 => undef,
    _seq2 => undef,
    _waiting => 1,
    _cmd => $cmdline
  };
  
  bless($self, $class);

  return $self;
}

sub destructor
{
  my ($self) = @_;

  close($self->{_in});
  close($self->{_out});

  waitpid($self->{_pid}, 1);
}

sub read_line
{
  my ($self) = @_;

  #print "SW Waiting...\n";

  my $in = $self->{_in};
  my $next_line;

  # Reading with time out
  # http://www.perlmonks.org/?node_id=43304
  eval {
    local $SIG{ALRM} = sub { die "timeout\n" };
    alarm($self->{_timeout});
    $next_line = <$in>;
  };

  if($@ eq "timeout\n") { die("Error: timeout reading NW output"); }
  elsif($@) { die("Error: couldn't read output"); }
  alarm(0);

  if(defined($next_line))
  {
    chomp($next_line);
    #print "IN: '$next_line'\n";

    if($next_line =~ /^Error/i)
    {
      print STDERR "ErrSeq1: '$self->{'seq1'}'\n";
      print STDERR "ErrSeq2: '$self->{'seq2'}'\n";
      croak($next_line);
    }
  }

  return $next_line;
}

sub read_line_fatal
{
  my ($self,$pattern) = @_;
  my $line = $self->read_line();
  if(!defined($line)) {
    croak("Error: cannot read from aligner.  Have you compiled?");
  } elsif($line !~ /$pattern/) {
    croak("Error: aligner said: $line (expected: $pattern)");
  }
}

sub do_alignment
{
  my ($self, $seq1, $seq2) = @_;

  my $line;
  my $out = $self->{_out};

  # _waiting is true if we have see the '==' line, indicating that the C
  # program is waiting for our sequence input
  if(!$self->{_waiting})
  {
    # Skip hits from previous alignment
    print $out "a\n";
    $line = $self->read_line_fatal(quotemeta(PROMPT_LINE."=="));
    $self->{_waiting} = 1;
  }

  if(length($seq1) == 0 || length($seq2) == 0)
  {
    carp("Cannot align lengths of zero");
    return;
  }
  elsif($seq1 =~ /[\n\r]/ || $seq2 =~ /[\n\r]/)
  {
    print STDERR "ErrSeq1: '$self->{'seq1'}'\n";
    print STDERR "ErrSeq2: '$self->{'seq2'}'\n";
    croak("New lines not allowed in sequences");
  }

  $self->{_align_number}++;
  $self->{_seq1} = $seq1;
  $self->{_seq2} = $seq2;

  my $expected = $self->{_align_number};

  print $out $seq1."\n";
  print $out $seq2."\n";

  $self->{_waiting} = 0;

  $line = $self->read_line_fatal("^== Alignment $expected");
  $line = $self->read_line_fatal('^$');
}

sub get_next_hit
{
  my ($self) = @_;

  if($self->{_waiting})
  {
    return undef;
  }

  # Print an 'h' to request next hit
  my $out = $self->{'_out'};
  print $out "h\n";

  my %result = ('seq1' => $self->{_seq1},
                'seq2' => $self->{_seq2});

  my $line;

  if(!defined($line = $self->read_line()))
  {
    die("No lines read in");
  }

  $line = substr($line, length(PROMPT_LINE));

  if($line =~ /^==/i)
  {
    # End of hits
    $self->{_waiting} = 1;
    return undef;
  }
  elsif($line =~ /^hit \d+\.(\d+) score: (\d+)$/i)
  {
    $result{'hit'} = $1;
    $result{'score'} = $2;

    my ($align1, $sep, $align2);

    if(defined($align1 = $self->read_line()) &&
       $align1 =~ /^  (.*)  \[pos: (\d+); len: (\d+)\]$/i)
    {
      $result{'align1'} = $1;
      $result{'pos1'} = $2;
      $result{'len1'} = $3;
    }
    else
    {
      die("Wasn't expecting '$align1'");
    }

    if(defined($sep = $self->read_line()) &&
       $sep =~ /^  ([\|\* ]+)$/)
    {
      $result{'sep'} = $1;
    }
    else
    {
      die("Wasn't expecting line '$sep'");
    }

    if(defined($align2 = $self->read_line()) &&
       $align2 =~ /^  (.*)  \[pos: (\d+); len: (\d+)\]$/i)
    {
      $result{'align2'} = $1;
      $result{'pos2'} = $2;
      $result{'len2'} = $3;
    }
    else
    {
      die("Wasn't expecting '$align2'");
    }

    # Skip remaining line
    $line = $self->read_line();

    return \%result;
  }
  else
  {
    die("Wasn't expecting '$line'");
  }
}

sub print_hit
{
  my ($self, $hit, $out) = @_;

  if(!defined($out))
  {
    open($out, ">-");
  }

  print $out "hit " .$self->{_align_number}.".".$hit->{'hit'}. " " .
             "score: " .$hit->{'score'}. ":\n";
  print $out "  " .$hit->{'align1'}. "  " .
             "[pos: " .$hit->{'pos1'}. "; len: " .$hit->{'len1'}. "]\n";
  print $out "  " .$hit->{'align2'}. "  " .
             "[pos: " .$hit->{'pos2'}. "; len: " .$hit->{'len2'}. "]\n";
}

1;
