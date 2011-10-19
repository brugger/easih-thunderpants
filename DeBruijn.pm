package DeBruijn;
# 
# General catch all and configuration wrapper module for all EASIH scripts.
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my %graph;

my $kmer = 27;
my $build_reads = 0;

my $EXIT_COUNTER = 20;

my $ID    = 0;
my $READS = 1;
my $IN    = 2;
my $OUT   = 3;

# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub set_kmer {
  my ($new_kmer) = @_;

  if ($new_kmer !~ /^\d+\z/ || 
      $new_kmer <= 0 ||
      $new_kmer >  35 ) {
    
    print STDERR "kmers should be a integer between 0 and 35 not '$new_kmer'\n";
    return;
  }


  $kmer = $new_kmer;
}



# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub dump_graph {
  print Dumper(\%graph);
}




# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub merge_singletons {

  my $merged_nodes = 0;

#  dump_graph();
# Remove edges only supported by a few reads
  my %handled = ('S' => 1);
 
  foreach my $node ( keys %graph ) {
    
    next if ( $handled{ $node } );

    # there is only one arch out from this node.
#    if ( keys %{$graph{ $node }{'OUT'} } == 1) {
    if ( _out( $node ) == 1) {
      my $edge  = (_out( $node ))[0];
      # The next node has more than one incoming arch
      next if ( $handled{ $edge });

      next if ( _in( $edge ) > 1);
      next if ( $handled{ $edge });

      _merge_nodes( $node, $edge);

      $handled{ $node }++;
      $handled{ $edge }++;
      $merged_nodes++;
    }
  }

#  print Dumper( \%graph );
  
  return $merged_nodes;
} 






# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub count_nodes {
  
  return keys %graph;
}




# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub _id {
  my ( $node, $id ) = @_;

  if ( $id ) {
    $graph{ $node }[ $ID ] = $id;
  }
       
  return $graph{ $node }[ $ID ];
}


# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub _reads {
  my ( $node, $reads ) = @_;

  if ( $reads && ref($reads) eq "HASH" ) {
    $graph{ $node }[ $READS ] = $reads;
  }
       
  return  \%{$graph{ $node }[ $READS ]};
}


# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub _in {
  my ( $node, $in ) = @_;

  if ( $in && ref($in) eq "HASH" ) {
    $graph{ $node }[ $IN ] = $in;
  }

  return keys %{$graph{ $node }[ $IN ]};
}


# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub _out {
  my ( $node, $out ) = @_;

  return keys %{$graph{ $node }[ $OUT]};
}



# 
# 
# 
# Kim Brugger (18 Oct 2011)
sub _merge_nodes {
  my ($node1, $node2) = @_;

  my $verbose_function = 0;

  if ( ! $graph{$node1} ) {
    print STDERR "$node1 does not exist cannot merge it with $node2\n";
    return;
  }

  if ( ! $graph{$node2} ) {
    print STDERR "$node2 does not exist cannot merge it with $node1\n";
    return;
  }

  my $merged_node = "";
  if ( substr($node2, -$kmer+1) eq substr($node1, 0, $kmer-1)) {
    ($node1, $node2) = ($node2, $node1);
  }

  if ( substr($node1, -$kmer+1) ne substr($node2, 0, $kmer-1)) {
    print STDERR "The two nodes does not overlap\n$node1  $node2\n";
    print STDERR substr($node1, -$kmer) ." ne ". substr($node2, 0, $kmer) ."\n";
    return;
  }

  # Calculate the new key
  $merged_node = $node1 . substr($node2, $kmer - 1 );    

  if ( _in( $node2 ) > 1 ) {
    print STDERR "$node2 has more than 1 incoming arc, cannot merge the node with $node1\n";
    return;
  }

  if ( _out( $node1 ) > 1 ) {
    print STDERR "$node1 has more than 1 outgoing arc, cannot merge the node with $node2\n";
    return;
  }

#  print STDERR  "New node: $merged_node\n";

  if ( $graph{ $merged_node } ) {
    print STDERR "$merged_node already exists (HOW DID THIS HAPPEN!!!)\n";
    return;
  }

  if ( $verbose_function ) {
    print "\n$node1  +  $node2 -->  $merged_node\n\n";
    print_node( $node1 );
    print_node( $node2 );
    print "\n";
  }

  # create the new node with the new key, and update the id
  $graph{ $merged_node } = $graph{ $node2 };

  _id( $merged_node, $merged_node);

  # Merge the counts of on the arches.
  my $reads1 = _reads( $node1 );
  my $reads2 = _reads( $merged_node );

  foreach my $read ( keys %$reads1 ) {
#    print STDERR Dumper($graph{ $node1 }{'OUT'}{$node2}{$read});
#    print STDERR Dumper($graph{ $node1 }{'OUT'}{$read});
    $$reads2{$read} ||= 0;
#    print STDERR " $read   $graph{ $merged_node }{'READS'}{$read} += $graph{ $node1 }{'READS'}{$read};\n";
    $$reads2{$read} += $$reads1{ $read };
  }

  _reads($merged_node, $reads2);

  # Move the incoming connections pointing to node1 onto the merged node
  _in( $merged_node, $graph{$node1}[ $IN ]);

  print "---- PRE [$merged_node] ----\n" if ( $verbose_function );
      
  # rename all the outgoing connections from the upstream node
  foreach my $pre ( _in( $merged_node )) {

    if ( ! $graph{ $pre }) {
      print STDERR "PRE :: $pre does not exist!\n";
      next;
    }
    
    print_node( $pre ) if ( $verbose_function );
    $graph{ $pre }[ $OUT ]{ $merged_node } = $graph{ $pre }[ $OUT ]{ $node1 };
    delete($graph{ $pre }[ $OUT ]{ $node1 });
    print_node( $pre ) if ( $verbose_function );
    
  }
  
  print "\n---- POST  [$merged_node ]----\n" if ( $verbose_function );
  # rename all the incoming connections from the upstream node
  foreach my $post ( _out( $merged_node )) {
    
    if ( ! $graph{$post}) {
      print STDERR  "POST :: $post does not exist!\n";
      next;
    }
    
    print_node( $post ) if ( $verbose_function );
    $graph{ $post }[ $IN ]{ $merged_node } = $graph{ $post }[ $IN ]{ $node2 };
    delete($graph{ $post }[ $IN ]{ $node2 });
    print_node( $post ) if ( $verbose_function );
    
  }
  
  delete( $graph{ $node1 });
  delete( $graph{ $node2 });      

#      print "----\n";
  if ( $verbose_function ) {
    print "\n";
    print_node( $merged_node );
  }
#  print STDERR Dumper ($graph{$merged_node});
#      print "-------\n\n";
      

}


# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub delete_low_weight {
  my ( $cutoff ) = @_;
  # Remove edges only supported by a few reads

  $cutoff ||= 0.1*$build_reads;


  my $purged = 0;
  foreach my $node ( keys %graph ) {
    my $reads = keys %{_reads( $node )};
    if ( $reads <= $cutoff ) {
      _delete_node( $node );
      $purged++;
    }
  }

  print STDERR "Removed $purged nodes with a $cutoff (total reads: $build_reads)\n";
  drop_orphans();

}



# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub print_tab {
  
  foreach my $node ( keys %graph ) {  
    foreach my $edge ( _out( $node )) {
      print join("\t", $node, $edge, $graph{$node}[ $OUT ]{$edge}, "\n");
    }
  }
}



# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub print_node {
  my ($key) = @_;
#  return;
  if ( ! _id( $key )) {
    print  "$key does not exist\n";
    return;
  }

  if ( ! _in( $key)) {
    print  "$key does not have any incoming arc\n";
#    return;
  }
  if ( ! _out( $key )) {
    print  "$key does not have any outgoing arc\n";
#    return;
  }

  print  join(" ", _in( $key )) . " --> $key --> " . join(" ", _out( $key )) . "\n";


}





# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub count_kmers {
  my ($filename) = @_;

#  print "Start counting\n";

  my %kmer_counts;
  my $reads;
  my ($seq);
  my $exit_counter = $EXIT_COUNTER;
  open( my $in, $filename) || die "Could not open '$filename': $!\n";
  while(<$in>) {
    chomp;

    if ( /\>/  ) {
      if ( $seq && length($seq) > 90) {
	for( my $i = 0; $i<length($seq) - $kmer; $i++) {
	  my $kmer = substr( $seq, $i, $kmer);
	  $kmer_counts{ $kmer }++;
	}
	
	$seq = "";
	$reads++;
	last if ($exit_counter-- == 0);
      }
    }
    else {
      $seq .= $_;
    }
  }
  for( my $i = 0; $i < length($seq) - $kmer; $i++) {
    my $kmer = substr( $seq, $i, $kmer);
    $kmer_counts{ $kmer }++;
  }
  close($in);
  $reads++;

#  print "Done counting\n";

  return ($reads, \%kmer_counts);
}



# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub readin_file {
  my ($filename) = @_;
  
  my ($reads, $kmer_counts);# = count_kmers( $filename );

  my ($name, $seq) = (1);
  my $exit_counter = $EXIT_COUNTER;

  open( my $in, $filename) || die "Could not open '$filename': $!\n";
  while(<$in>) {
    chomp;

    if ( /\>/  ) {
      if ( $seq && length($seq) > 90) {
#      print "$name\n";
	DeBruijn::add_sequence( $name, $seq, $kmer_counts, $reads );
	last if ( $exit_counter-- == 0 );
 	$seq = "";
      }
      $name++;
    }
    else {
      $seq .= $_;
    }
  }

  DeBruijn::add_sequence( $name, $seq, $kmer_counts, $reads ) if ( $seq );
#  print Dumper( \%graph);
}

# 
# 
# 
# Kim Brugger (06 Oct 2011)
sub add_sequence {
  my ($name, $seq, $kmer_counts, $reads) = @_;

#  print "$seq\n";

  $build_reads++;
  my $prev_node = "S";

  for( my $i = 0; $i<length($seq) - $kmer; $i++) {
    my $new_node = substr( $seq, $i, $kmer);


    if ( $kmer_counts && $reads && 
	 $$kmer_counts{ $new_node } < $reads*0.1 ){
      $prev_node = undef;
      next;
    }
	  
    if ( $prev_node ) {
#      print "$prev_node $new_node --> $$kmer_counts{ $new_node } $reads\n";

      $graph{ $prev_node }[ $READS ]{ $name }++;
      $graph{ $prev_node }[ $ID    ] = $prev_node;

      $graph{ $new_node  }[ $READS ]{ $name }++;
      $graph{ $new_node  }[ $ID    ] = $kmer;

      $graph{ $new_node  }[ $IN    ]{ $prev_node}++;
      $graph{ $prev_node }[ $OUT   ]{ $new_node }++;
    }
    $prev_node = $new_node;

  }
}


# 
# 
# 
# Kim Brugger (18 Oct 2011)
sub find_accessible_nodes {
  my ($node) = @_;

  $node ||= 'S';

  my @nodes =( $node );
  
  foreach my $out_node ( _out( $node ) ) {
    push @nodes, find_accessible_nodes( $out_node );
  }

  return @nodes;

}


# 
# 
# 
# Kim Brugger (18 Oct 2011)
sub drop_orphans {

  my %nodes = map{ $_ => 1 } keys %graph;

  my @live_nodes = find_accessible_nodes();

#  print Dumper( \@live_nodes );
#  print Dumper( \%nodes );

  foreach my $node ( @live_nodes ) {
    delete($nodes{ $node});
  }

  print STDERR "Dropping " . (keys %nodes) . " orphans \n";

  foreach my $node ( keys %nodes ) {
    _delete_node($node);
  }
  

}



# 
# 
# 
# Kim Brugger (18 Oct 2011)
sub _delete_node {
  my ( $node ) = @_;

  foreach my $in ( _in( $node ) ) {
    delete ($graph{$in}[$OUT]{ $node });
  }

  foreach my $out ( _out( $node )) {
    delete ($graph{ $out }[ $IN ]{ $node });
  }

  delete( $graph{$node} );

  
}



# 
# 
# 
# Kim Brugger (07 Oct 2011)
sub path_finder {


  my $pos = 'S';


#  print "start paths: ". join(" ", keys %{$graph{'S'}{'OUT'}}) . "\n";

#  print Dumper(\%graph);

  foreach my $start_pos ( keys %{$graph{'S'}{'OUT'}} ) {
#    print "S $start_pos\n";
    my @post_poss = keys %{$graph{$start_pos}{'OUT'}};
    foreach my $post_pos ( @post_poss ) {

#      print "S $start_pos $post_pos\n";
      my $shared_read = 0;
      foreach my $read ( keys %{$graph{ $start_pos }{'OUT'}{$post_pos}} ) {
	if ( $graph{ $pos }{'OUT'}{$start_pos}{$read}) {
	  print "SHARED READ S:: $start_pos -- $post_pos $read :: \n";
#	  print "[$graph{ 'S' }{'OUT'}{$start_pos}{$read} --";
#	  print "$graph{ $start_pos }{'OUT'}{$post_pos}{$read}]\n";
	  $shared_read++;
	  last;
	}
      }

      if ($shared_read ) {
	my %new_legacy = (%{$graph{ 'S' }{'OUT'}{$start_pos}}, %{$graph{ $start_pos }{'OUT'}{$post_pos}});
	_path_finder($start_pos.substr($post_pos, $kmer -1), $post_pos, \%new_legacy);
      }
      else {
	print "Dropping path $start_pos $post_pos\n";
      }
    }
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


  my @post_poss = keys %{$graph{$pos}{'OUT'}};
#  if (! @post_poss  || @post_poss > 1) {
  if (! @post_poss ) {
    my $weight = keys %$legacy;
#    print "PATH :: $pre_path\t$weight\n" if ( length($pre_path) > 90);
    print "PATH :: $pre_path\t$weight\n" ;
    return;
  }

  foreach my $post_pos ( @post_poss ) {

    print "Trying to connect with $post_pos\n";
    my $shared_read = 0;
    foreach my $read ( keys %{$graph{ $pos }{'OUT'}{$post_pos}} ) {
      if ( $$legacy{$read}) {
	  print "SHARED READ :: $pos -- $post_pos $read :: \n";
#	  print "[$graph{ $pos }{'OUT'}{$post_pos}{$read} --";
#	  print "$$legacy{$read}]\n";
	$shared_read++;
	last;
      }
    }
    
    
    if ( $shared_read ) {
#      my %new_legacy = (%{$graph{ 'S' }{'OUT'}{$start_pos}}, %{$graph{ $start_pos }{'OUT'}{$post_pos}});
      my %new_legacy = (%{$graph{ $pos }{'OUT'}{$post_pos}}, %$legacy);
      _path_finder ($pre_path .substr($post_pos, $kmer -1), $post_pos, \%new_legacy );
    }
    else {
      print "Dropping path $pos $post_pos\n";
    }
  }
  
  
}



1;



__END__


# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub merge_singletons_old {

  my $merged_nodes = 0;

# Remove edges only supported by a few reads
  my %handled = ('E' => 1);
 
  foreach my $node ( keys %graph ) {
    
    next if ( $handled{ $node } );

    # there is only one arch out from this node.
    if ( keys %{$graph{ $node }{'OUT'} } == 1) {
      my $edge  = (keys %{$graph{ $node }{'OUT'}})[0];
      # The next node has more than one incoming arch


      next if ( keys %{$graph{ $edge }{'IN'}} > 1);

      next if ( $handled{ $edge });

      # Calculate the new key
      my $merged_node = $node . substr($edge, $kmer - 1 );

      if ( $graph{ $merged_node } ) {
	print "$merged_node already exists\n";
	next;
      }

#      print "$node  +  $edge -->  $merged_node\n";
#      print_node( $node );
#      print_node( $edge );

      # create the new node with the new key
      $graph{ $merged_node } = $graph{ $edge };
      %{$graph{ $merged_node }{'IN'}} = %{$graph{ $node }{'IN'}} if ($graph{ $node }{'IN'});

#      print "---- PRE ----\n";
      
      # rename all the Incomming connections from the upstream node
      foreach my $pre ( keys %{$graph{ $merged_node }{'IN'} }) {

	if ( ! $graph{$pre}) {
	  print STDERR "PRE :: $pre does not exist!\n";
	  next;
	}
#	print_node( $pre );
	$graph{ $pre }{ 'OUT' }{ $merged_node } = $graph{ $pre }{ 'OUT' }{ $node };
	delete($graph{ $pre }{ 'OUT' }{ $node });
#	print_node( $pre );
	
      }

#      print "---- POST ----\n";
      # rename all the Incomming connections from the upstream node
      foreach my $post ( keys %{$graph{ $merged_node }{'OUT'} }) {

	if ( ! $graph{$post}) {
	  print  "POST :: $post does not exist!\n";
	  next;
	}

#	print_node( $post );
	$graph{ $post }{ 'IN' }{ $merged_node } = $graph{ $post }{ 'IN' }{ $edge };
	delete($graph{ $post }{ 'IN' }{ $edge });
#	print_node( $post );
	
      }
      
      delete( $graph{ $node });
      delete( $graph{ $edge });      

#      print "----\n";
#      print_node( $merged_node );
#      print "-------\n\n";
      

      $handled{ $node }++;
      $handled{ $edge }++;
      $merged_nodes++;
    }
  }

#  print Dumper( \%graph );
  
  return $merged_nodes;
} 




# 
# 
# 
# Kim Brugger (06 Oct 2011)
sub add_sequence_old {
  my ($name, $seq, $kmer_counts, $reads) = @_;

#  print "$seq\n";
  
  my $prev_node = "S";

  for( my $i = 0; $i<length($seq) - $kmer; $i++) {
    my $new_node = substr( $seq, $i, $kmer);


    if ( $kmer_counts && $reads && 
	 $$kmer_counts{ $new_node } < $reads*0.1 ){
      $prev_node = undef;
      next;
    }
	  
    if ( $prev_node ) {
#      print "$prev_node $new_node --> $$kmer_counts{ $new_node } $reads\n";

      $graph{ $new_node }{'READS'}{ $name }++;
      $graph{ $new_node }{'ID'   } = $kmer;

      $graph{ $prev_node }{'OUT' }{ $new_node  }{$name}++;
      $graph{ $new_node }{ 'IN'  }{ $prev_node }{$name}++;
    }
    $prev_node = $new_node;

  }
}


