package SmithWaterman;

use strict;
use warnings;

use constant DEFAULT_CMD => 'smith_waterman';

use Carp;
use FileHandle;

use IPC::Open2;
use IPC::Open3;
use IO::Select;

#use base 'Exporter';
#our @EXPORT = qw(method_name);

sub new
{
  my ($class, @args) = @_;

  my %options = (@args);

  my $cmdline = DEFAULT_CMD;

  for my $key (keys %options)
  {
    $options{lc($key)} = $options{$key};
  }

  if(defined($options{'mismatch'}) != defined($options{'mismatch'}))
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
    if($key =~ /^case_sensitive$/i)
    {
      if($value)
      {
        $cmdline .= " --case_sensitive";
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
    elsif(grep(/^$key$/i, qw(substitution_matrix substitution_pairs
                             match mismatch gapopen gapextend
                             minscore maxhits nogaps nomismatches wildcard)))
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
  #use Symbol 'gensym'; $err = gensym;
  #my $pid = open3($out, $in, $err, $cmdline)
  #  or die("Cannot run cmd: '$cmdline'");

  print "running '$cmdline'\n";

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
    _waiting => 1,
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

  my $in = $self->{_in};
  my $next_line = <$in>;

  chomp($next_line);
  #print "line: '$next_line'\n";

  #my $next_line;
  #my $in = $self->{_in};
  #my $err = $self->{_err};

  #my $select = new IO::Select();

  #$select->add($in);
  #$select->add($err);

  #my $timeout = 1;

  #foreach my $h ($select->can_read($timeout))
  #{
  #  if($h eq $err)
  #  {
      #sysread($err,$buf,4096);
     #if($buf){print "ERROR-> $buf\n"}
  #    print "Reading STDERR\n";
  #    $next_line = <$err>;
  #    croak("Error: $next_line\n");
  #  }
  #  else
  #  {
     #sysread($in,$buf,4096);
     #if($buf){print "in: $buf\n"}
  #    print "Reading STDOUT\n";
  #    $next_line = <$in>;
  #    chomp($next_line);
  #    print "STDOUT: '$next_line'\n";
  #  }
  #}

  if(defined($next_line) && $next_line =~ /^SmithWaterman Error/i)
  {
    croak($next_line);
  }

  return $next_line;
}

sub do_alignment
{
  my ($self, $seq1, $seq2) = @_;

  my $line;

  if(!$self->{_waiting})
  {
    # Skip hits from previous alignment
    while(defined($line = $self->read_line()) && $line =~ /^hit/i)
    {
      for(my $i = 0; $i < 3; $i++)
      {
        $line = $self->read_line();

        if($line !~ /^  /)
        {
          die("Unexpected output '$line'");
        }
      }

      $line = $self->read_line();

      if($line !~ /^$/)
      {
        die("Unexpected output '$line'");
      }
    }

    if($line ne "==")
    {
      die("Unexpected output '$line'");
    }

    $self->{_waiting} = 1;
  }

  if(length($seq1) == 0 || length($seq2) == 0)
  {
    carp("Cannot align lengths of zero");
    return;
  }
  elsif($seq1 =~ /[\n\r]/ || $seq2 =~ /[\n\r]/)
  {
    croak("New lines not allowed in sequences");
  }

  my $in = $self->{_in};
  my $out = $self->{_out};

  $self->{_align_number}++;
  $self->{_seq1} = $seq1;
  $self->{_seq2} = $seq2;

  my $expected = $self->{_align_number};

  print $out $seq1."\n";
  print $out $seq2."\n";

  $self->{_waiting} = 0;

  if(!defined($line = $self->read_line()) ||
     $line !~ /^== Alignment $expected/i)
  {
    die("Wasn't expecting '$line'");
  }

  if(!defined($line = $self->read_line()) ||
     $line !~ /^$/i)
  {
    die("Wasn't expecting '$line'");
  }
}

sub get_next_hit
{
  my ($self) = @_;

  if($self->{_waiting})
  {
    return undef;
  }

  my %result = ('seq1' => $self->{_seq1},
                'seq2' => $self->{_seq2});

  my $in = $self->{_in};
  my $line;

  if(!defined($line = $self->read_line()))
  {
    die("No lines read in");
  }
  elsif($line =~ /^==/i)
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

sub sw_print_hit
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
