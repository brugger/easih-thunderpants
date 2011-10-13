#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (06 Oct 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

no warnings "recursion";


my %graph;
my %starts;
my $kmer = 27;


my $seq;

my $exit_counter = 200;

my $reads = 0;
my $name = 1;

while(<>) {
  chomp;

  if ( /\>/  ) {
    if ( $seq ) {
#      print "$name\n";
      add_to_graph( $name, $seq );
      last if ($exit_counter-- <= 0);
      $seq = "";
      $reads++;
    }
    $name++;
#    $name = $_;
  }
  else {
    $seq .= $_;
  }
}


foreach my $sub_seq ( keys %graph ) {
  
  foreach my $edge (keys %{$graph{ $sub_seq } } ) {

    if ( keys %{$graph{ $sub_seq }{ $edge }} <= 20) {
#      print "Deleting $sub_seq --> $edge ".(keys %{$graph{ $sub_seq }{ $edge }})."\n";
#      delete( $graph{$edge});
      delete($graph{ $sub_seq }{ $edge });
    }
  }

  delete( $graph{ $sub_seq } ) if ( keys %{$graph{ $sub_seq } } == 0 );
  
}



print Dumper( \%starts );
print Dumper( \%graph );
#exit;

print "-------------\n";


#print Dumper(\%starts );

foreach my $start ( keys %starts ) {  

  #print Dumper( $graph{ $start}) if ($graph{ $start});

  next if ( $starts{$start} < 2 );

  path_finder("", $start );  
  next;
  my %paths;


  if ($graph{ $start}) {
    my $point = $start;
    print $start;
    my $weight = 0;
    while ( 1 ) {
      my @keys = keys %{$graph{$point}};
#      print "\n" . @keys . " $point";
      last if ( ! @keys ); 
#      @keys = sort { @{$graph{$point}{$a}} <=> @{$graph{$point}{$b}}} @keys;
#      print "\n-" . substr( $keys[0], -1) . "($keys[0]) " . (@keys) . " " ;
      print substr( $keys[-1], -1);
      $weight += @{ $graph{$point}{ $keys[-1] }};
      $point = $keys[-1];
    }

    print " $weight\n";
#    exit;

  }
}



# 
# 
# 
# Kim Brugger (07 Oct 2011)
sub path_finder {
  my ($pre_path, $pos) = @_;

  my @post_poss = keys %{$graph{$pos}};
  foreach my $post_pos ( @post_poss ) {
#    _path_finder($post_pos, $post_pos);
    _path_finder($pos. substr($post_pos, -1), $post_pos, $graph{ $pos }{$post_pos});
  }
  
  
}


# 
# 
# 
# Kim Brugger (07 Oct 2011)
sub _path_finder {
  my ($pre_path, $pos, $legacy) = @_;

#  die Dumper( $legacy );

#  print "$pre_path -- $pos\n";


  my @post_poss = keys %{$graph{$pos}};
#  if (! @post_poss  || @post_poss > 1) {
  if (! @post_poss ) {
    print "$pre_path\n";
    return;
  }
  foreach my $post_pos ( @post_poss ) {
    
    my $shared_read = 0;
    foreach my $read ( keys %{$graph{ $pos }{$post_pos}} ) {
      if ( $$legacy{ $read }) {
#      print "SHARED READ :: $read $post_pos\n";
	$shared_read++;
	last;
      }
    }

    
    

    next if ( !$shared_read );
    
    if ( $shared_read ) {
      my %new_legacy = (%{$graph{ $pos }{$post_pos}}, %$legacy);
      _path_finder ($pre_path ."" .substr($post_pos, -1), $post_pos, \%new_legacy );
    }
  }
  
  
}



#print Dumper( \%starts );

#print "Readin $reads reads\n";
#print Dumper( \%graph );


# 
# 
# 
# Kim Brugger (06 Oct 2011)
sub add_to_graph {
  my ($name, $seq) = @_;

#  print "$seq\n";
  
  my $old_kmer = "";

  for( my $i = 0; $i<length($seq) - $kmer; $i++) {
    my $sub_seq = substr( $seq, $i, $kmer);
#    print "$sub_seq $old_kmer\n";
    if ( $old_kmer) {
#      $graph{$old_kmer}{$sub_seq}++;
      $graph{$old_kmer}{$sub_seq}{$name}++;
      if ( $i == 1 ) {
	$starts{ $old_kmer }++;
      }
    }
    $old_kmer = $sub_seq;
  }

}
