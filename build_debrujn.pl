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

use lib './';
use lib '../';
use lib '../../';
use DeBruijn;


my %graph;
my $kmer = 27;


my $seq;

my $exit_counter = 200;

my $reads = 0;
my $name = 1;

if (0) {
  while(<>) {
    chomp;
    
    if ( /\>/  ) {
      if ( $seq ) {
#      print "$name\n";
	DeBruijn::add_sequence( $name, $seq );
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
}
else {
  DeBruijn::readin_file(shift);
}

#DeBruijn::_merge_nodes('CACATTTCTTGGAGTACTCTACGTCTGAGT', 'TTTCTTGGAGTACTCTACGTCTGAGTG');
#DeBruijn::_merge_nodes('CAATGGGACGGAGCGGGTGCGGTTCCT', 'AATGGGACGGAGCGGGTGCGGTTCCTG');
#exit;

DeBruijn::delete_low_weight();
DeBruijn::drop_orphans();



#DeBruijn::delete_low_weight();
#DeBruijn::merge_singletons();
print STDERR "Start simplifying graph with ".DeBruijn::count_nodes()." nodes...\n";
while (DeBruijn::merge_singletons() ) {
  print STDERR "Simplifying graph (".DeBruijn::count_nodes()." nodes)...\n";
  ;
}
print STDERR "Done simplifying graph, ".DeBruijn::count_nodes()." nodes left\n";
DeBruijn::drop_orphans();



#DeBruijn::_merge_nodes('TTTCTTGGAGTACTCTACGTCTGAGTG', 'CACATTTCTTGGAGTACTCTACGTCTGAGT');


DeBruijn::print_tab();
DeBruijn::path_finder();
#DeBruijn::dump_graph();



foreach my $start ( keys %graph ) {  

  foreach my $key (keys %{$graph{ $start }} ) {

    print join("\t", $start, $key, int(keys %{$graph{$start}{$key}}), "\n");

  }
}

exit;



#print Dumper(\%starts );

foreach my $start ( keys %{$graph{S}} ) {  

  #print Dumper( $graph{ $start}) if ($graph{ $start});


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


