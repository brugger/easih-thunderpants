#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (04 Nov 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


use lib './';
use lib '../';
use lib '../../';
use Typer;

my $file = shift;
Typer::readin_file($file) if ($file =~ /.fa/);
Typer::readin_db($file)   if ($file =~ /.tkr/);

#my @hits;
my %hits;
my $fas = readin_fasta_file( shift );

#print Dumper( $fas );

my $region_index = 0;

my %type_index;
foreach my $region ( sort keys %$fas ) {

  print "\n\nRegion == $region\n";
  
  my (@tmp_scores, $best_score);
  foreach my $seq ( @{$$fas{ $region}} ) {
    
    my ($score, $types) = Typer::type_sequence( $seq );
#    @{$types} = grep(/0[37]:01/, @{$types});
#    print Dumper( $types );
    my $type_head = substr(join(" ", @{$types}), 0 , 50);
    $type_head = join(" ", @{$types});
    print "SCORE :: $score, $type_head\n";
#    next if ( $score < 80);
#    next if ( $type_head !~ /DRB1/);
#    next if ( $type_head =~ /DRB3/);
    push @tmp_scores, [$score, $types];
    $best_score = $score if ( ! $best_score || $best_score < $score );
  }
  
  @tmp_scores = sort { $$b[0] <=> $$a[0] } @tmp_scores;

#  print Dumper( \@tmp_scores );

  my $picked = 0;
  my $last_score;
  my $read_index = 0;
  foreach my $score ( @tmp_scores ) {
    
    $last_score = $$score[ 0 ] if ( !$last_score );
    last if ( $last_score != $$score[ 0 ] && $picked >= 2);
    
#    print Dumper( $score );

    foreach my $type ( @{$$score[1]}) {
#      push @{$hits{ $type }{ $region  }}, $read_index;
      $hits{ $type }{ $region  }{$read_index}++;
    }
    $picked++;

    $last_score = $$score[ 0 ];
    $read_index++;
  }


  $region_index++;
}   

my $exit_count = 30;
foreach my $coverage ( sort { keys %{$hits{ $b }} <=> keys %{$hits{ $a }}} keys %hits ) {
  print "$coverage :: ". (keys %{$hits{ $coverage }}). "\n";
  last if ( $exit_count-- <= 0);
}

print "\n\n\n";
#print Dumper( \%hits );

my $printed = 0;

 REDO:
my %support_reads = ();

my $score = 0;
foreach my $coverage ( sort { keys %{$hits{ $b }} <=> keys %{$hits{ $a }}} keys %hits ) {
  $score = keys %{$hits{ $coverage }} if ( ! $score );
  last if (  $score != keys %{$hits{ $coverage }} );

  print "$coverage :: ". (keys %{$hits{ $coverage }}). " $printed\n";
  foreach my $region (keys %{$hits{ $coverage }}) {
    foreach my $read_index ( keys %{$hits{ $coverage }{$region}} ){
      $support_reads{ $region }{ $read_index }++;
    }
  }

  
}
$printed++;

if ( $printed < 2 ) {
  foreach my $type (keys %hits) {
    foreach my $region (keys %{$hits{ $type }}) {
      foreach my $read_index ( keys %{$hits{ $type }{$region}} ){
	delete $hits{ $type }{$region}{ $read_index } if ( $support_reads{$region}{ $read_index} );
      }
      delete $hits{ $type }{$region} if ( ! keys %{$hits{ $type }{$region}});
    }
    delete $hits{ $type } if ( ! keys %{$hits{ $type }});
  }
  
  goto REDO ;
}

print "\n\n";
#print Dumper( \%hits );

#foreach my $type ( keys %hits  ) {
  
#  delete $hits{$type} if ( keys %{$hits{$type}} < 8 );

#}

#print Dumper( \%hits );

# 
# 
# 
# Kim Brugger (06 Nov 2011)
sub readin_fasta_file {
  my ( $filename ) = @_;  

  my %res;
  my ($name, $seq) = ("", "");
  open( my $in, $filename) || die "Could not open '$filename': $!\n";
  while(<$in>) {
    chomp;
    
    if ( /\>/  ) {
      if ( $seq ) {
	push @{$res{$name}}, $seq;
      }
      $seq = "";
      $name = $_;
      $name =~ s/\>//;
      $name =~ s/.fa//;
#      $name =~ s/_[FR]//;

    }
    else {
      $seq .= $_;
    }
  }
  if ( $name && $seq ) {
    push @{$res{$name}}, $seq;
  }
  
  return \%res;
}

