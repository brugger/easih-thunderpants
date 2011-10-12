#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (11 Oct 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $ref = shift;
$ref = readin_ref( $ref );

my %stats;

my (%hist);

while(<>) {
  chomp;
  my @F = split("\t");
  my ($name, $flags, $pos, $cigar, $sequence, $qual) = ($F[0], $F[1], $F[3], $F[5], $F[9], $F[10]);

  my ($psequence, $pqual) = patch_alignment($sequence, $qual, $cigar);


  my $mask = alignment_mask($cigar);
  my $sub_ref = substr( $ref, $pos - 1, length( $psequence));

#  print "$sub_ref\n$psequence\n$pqual\n";

  my @sub_ref = split(//, $sub_ref);
  my @psequence = split(//, $psequence);
  my @pqual = split(//, $pqual);

  for(my $i = 0; $i< @sub_ref; $i++) {

    $hist{ord($pqual[$i])-33}{M}++;

    next if ( $psequence[ $i ] eq "-" ||
	      $sub_ref[$i]     eq "-");    

    if ( $sub_ref[$i] eq $psequence[ $i ] ) {
      $stats{ ord($pqual[$i]) - 33}{$psequence[ $i ]}{ M }++;
      $stats{ ord($pqual[$i]) - 33}{ M }++;
    }
    else {
      $stats{ ord($pqual[$i]) - 33 }{$psequence[ $i ]}{ X }++;
      $stats{ ord($pqual[$i]) - 33 }{ X }++;
    }

  }

  foreach my $QV ( split(//, $qual )) {
    $hist{ ord($QV) - 33}{R}++;
  }
  
}

foreach my $QV ( sort {$a <=> $b} keys %stats ) {

  
  printf("$QV\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\n",
	 phred_value( $stats{ $QV }{ M }, $stats{ $QV }{ X }),
	 phred_value( $stats{ $QV }{ A }{ M }, $stats{ $QV }{ A }{ X }),
	 phred_value( $stats{ $QV }{ C }{ M }, $stats{ $QV }{ C }{ X }),
	 phred_value( $stats{ $QV }{ G }{ M }, $stats{ $QV }{ G }{ X }),
	 phred_value( $stats{ $QV }{ T }{ M }, $stats{ $QV }{ T }{ X }),
	 $hist{$QV}{M} || 0,
	 $hist{$QV}{R} || 0,

);



}



# 
# 
# 
# Kim Brugger (11 Oct 2011)
sub phred_value {
  my ($correct, $wrong) = @_;
  
  return -1 if ( !$wrong && ! $correct);
  return 55 if ( ! $wrong && $correct);
  return  0 if ( $wrong && ! $correct);

  $wrong ||= 0;

  my $P = $wrong/($wrong+$correct);
  
  return 55 if ( ! $P );
  my $Q = -10 * log10( $P );
  return int($Q);
}




sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

# 
# 
# 
# Kim Brugger (10 Oct 2011)
sub explode_sequence {
  my ($seq) = @_;
  
  my @stack;
  my $last;
  my $bla = "";
  foreach my $base ( split(//, $seq) ) {
    if ($last && $last ne $base ) {
      push @stack, $bla;
      $bla = "";
    }
    $bla .= $base;
    $last = $base;
  }

  push @stack, $bla;
  
#  print Dumper( \@stack );
  return \@stack;

  exit;

}



# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub alignment_mask {
  my ($cigar ) = @_;

  my $mask = "";

  my (@cigar) = $cigar =~ /(\d*\w)/g;


  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)



  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d*)(\w)/;
    $length ||= 1;
    
    $mask .= "$type"x$length;
  }
    return $mask;
}



# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub patch_alignment {
  my ( $seq, $qual, $cigar ) = @_;

  return ($seq, $qual) if ( $cigar !~ /[DIS]/);
  
  my @seq  = split("", $seq );
  my @qual = split("", $qual );


  my (@cigar) = $cigar =~ /(\d*\w)/g;

  my $offset = 0;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)

  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d+)(\w)/;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "D") {
      my @dashes = split("", "-"x$length);
      splice(@seq,  $offset, 0, @dashes);
      splice(@qual, $offset, 0, @dashes);
      $offset += $length;
    }
    elsif ( $type eq "I" || $type eq "S" ) {
      splice(@seq,  $offset, $length);
      splice(@qual, $offset, $length);
    }    

  }


  return (join("", @seq), join("",@qual));
}


# 
# 
# 
# Kim Brugger (20 Jul 2009)
sub patch_alignment2 {
  my ( $seq, $qual, $cigar ) = @_;

  return ($seq, $qual) if ( $cigar !~ /[DIS]/);
  
  my @seq  = split("", $seq );
  my @qual = split("", $qual );


  my (@cigar) = $cigar =~ /(\d*\w)/g;

  my $offset = 0;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)



  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d+)(\w)/;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "D") {
      my @dashes = split("", "-"x$length);
      splice(@seq,  $offset, 0, @dashes);
      splice(@qual, $offset, 0, @dashes);
      $offset += $length;
    }
    elsif ( $type eq "I" || $type eq "S" ) {
      splice(@seq,  $offset, $length);
      splice(@qual, $offset, $length);
    }    

  }


  return (join("", @seq), join("",@qual));
}



# 
# 
# 
# Kim Brugger (10 Oct 2011)
sub readin_ref {
  my ($file) = @_;

  my $seq;
  open(my $in, $ref) || die "Could not open '$file':$!\n";
  while(<$in>) {
    next if ( /^\>/);
    chomp;
    $seq .= $_;
  }
  close( $in );
  return $seq;
}
