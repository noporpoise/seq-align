package NeedlemanWunsch;

use strict;
use warnings;

use constant DEFAULT_CMD => 'needleman_wunsch';

use Carp;
use FileHandle;

use IPC::Open2;

sub new
{
  my ($class, @args) = @_;

  my %options = (@args);

  my $cmdline = DEFAULT_CMD;

  for my $key (keys %options)
  {
    $options{lc($key)} = $options{$key};
  }

  if(defined($options{'match'}) != defined($options{'mismatch'}))
  {
    carp("Cannot set only one of match/mismatch");
  }

  if(defined($options{'cmd'}))
  {
    $cmdline = $options{'cmd'};
  }

  $cmdline .= " --stdin --pretty --printscores";

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
        carp("Invalid NeedlemanWunsch.pm scoring '$value'");
      }
      else
      {
        $cmdline .= " --scoring $value";
      }
    }
    elsif(grep(/^$key$/i, qw(substitution_matrix substitution_pairs
                             match mismatch gapopen gapextend wildcard)))
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

  $out->autoflush();

  my $self = {
    _in => $in,
    _out => $out,
    _err => $err,
    _pid => $pid,
    _align_number => -1,
    _seq1 => undef,
    _seq2 => undef,
    _cmd => $cmdline
  };
  
  bless $self, $class;

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
  my $next_line = <$in>;

  chomp($next_line);
  #print "IN: '$next_line'\n";

  if(defined($next_line) && $next_line =~ /^NeedlemanWunsch Error/i)
  {
    croak($next_line);
  }

  return $next_line;
}

sub do_alignment
{
  my ($self, $seq1, $seq2) = @_;

  if($seq1 =~ /[\n\r]/ || $seq2 =~ /[\n\r]/)
  {
    croak("New lines not allowed in sequences");
  }

  my $out = $self->{_out};

  print $out "$seq1\n";
  print $out "$seq2\n";

  my %result = ('seq1' => $seq1, 'seq2' => $seq2,
                'number' => $self->{_align_number}++);

  $result{'align1'} = $self->read_line();
  $result{'sep'} = $self->read_line();
  $result{'align2'} = $self->read_line();
  my $score_line = $self->read_line();

  if(!defined($result{'align1'} || !defined($result{'sep'}) ||
     !defined($result{'align2'}) || !defined($score_line)))
  {
    die("Missing lines when reading in");
  }

  if($score_line =~ /(-?\d+)$/i)
  {
    $result{'score'} = $1;
  }
  else
  {
    croak("Cannot locate score in string '".$score_line."'");
  }

  # Skip remaining line
  my $skip = $self->read_line();
  chomp($skip);

  if($skip ne '')
  {
    croak("Not expecting line: '$skip'");
  }

  return \%result;
}

sub print_alignment
{
  my ($self, $hit, $out) = @_;

  if(!defined($out))
  {
    open($out, ">-");
  }

  print $out $hit->{'align1'}."\n";
  print $out $hit->{'sep'}."\n";
  print $out $hit->{'align2'}."\n";
  print $out "score: $hit->{'score'}\n\n";
}

1;
