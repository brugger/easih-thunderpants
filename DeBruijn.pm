package DeBruijn;
# 
# General catch all and configuration wrapper module for all EASIH scripts.
# 
# 
# Kim Brugger (13 Jun 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
no warnings "recursion";

my %graph;

my $kmer = 12;
my $build_reads = 0;

my $MIN_RATIO = 0.10;

my $EXIT_COUNTER = -4000;

my $ID    = 0;
my $READS = 1;
my $IN    = 2;
my $OUT   = 3;

my $REVERSE = 0;

my $FRAG_MIN = 90;
my $FRAG_MAX = 129;


my $infile = "";

my %multipliers; # identical sequences get stored 

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
# Kim Brugger (28 Oct 2011)
sub max_length {
  my ( $node, $length) = @_;
  $length ||= $kmer;
  $length += length( $node ) - $kmer;

  my @starts = _out($node);
  my $max_length = $length;

  foreach my $start ( @starts ) {
    my $cur_length = max_length($start, $length);
    $max_length = $cur_length if ( $max_length < $cur_length );
  }
  
  return $max_length;
}





# 
# 
# 
# Kim Brugger (28 Oct 2011)
sub collapse_start_frags {

  my @starts = _out('S');
  my %starts_hash = map{ $_ => 1} @starts;
  foreach my $start ( @starts ) {
    my $max_length = max_length( $start );
    
    if ( $max_length < $FRAG_MIN ) {
      print "$start will generate a fragment of max: " . max_length( $start ) . " bp, deleting it\n";
      _delete_node( $start );
    }

    my @downstreamed =  _downstream_members($start, 2);
    foreach my $down ( @downstreamed ) {
      next if ( $start eq $down );
      if ( $starts_hash{ $down } ) {
	print "$start could be joined with $down\n";
      }
    }
  }
  drop_orphans();
#  print join("\n", @starts, "\n");

}





# 
# 
# 
# Kim Brugger (28 Oct 2011)
sub _downstream_members {
  my ($node, $depth) = @_;
  return if ($depth <= 0);

  my @kids = ($node);
  foreach my $out (_out($node)) {
    push @kids, _downstream_members($out, $depth - 1);
  } 

  return @kids;
}



# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub merge_singletons {

  my $merged_nodes = 0;

  my %handled = ('S' => 1);
 
  foreach my $node ( keys %graph ) {
    
    next if ( $handled{ $node } );

    # there is only one arch out from this node.
#    if ( keys %{$graph{ $node }{'OUT'} } == 1) {
    if ( _out( $node ) == 1) {
      my $edge  = (_out( $node ))[0];
      next if ( $handled{ $edge });

      # The next node has more than two incoming arch
      next if ( _in( $edge ) > 2);
      #, check and see if one is a shared 'S' node, if so merge the bastards!
      if ( _in( $edge ) == 2) {
	my $node_in = grep(/^S\z/, _in( $node ));
	my $edge_in = grep(/^S\z/, _in( $edge ));
	
	next if (  ! $node_in || ! $edge_in );
	delete($graph{'S'}[ $OUT ]{ $edge });
	delete($graph{$edge}[ $IN ]{ 'S' });

      }

	
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

  $cutoff ||= $MIN_RATIO*$build_reads;

  print STDERR  " ".DeBruijn::count_nodes()." nodes before pruning...";

  my $purged = 0;
  foreach my $node ( keys %graph ) {
    my $reads = keys %{_reads( $node )};
    if ( $reads <= $cutoff ) {
      _delete_node( $node );
      $purged++;
    }
  }
  print STDERR   " ".DeBruijn::count_nodes()." nodes left\n";

  print STDERR "Removed $purged nodes with a $cutoff depth cutoff (total reads: $build_reads)\n";
  drop_orphans();
}



# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub print_tab {
  
  foreach my $node ( keys %graph ) {  
    foreach my $edge ( _out( $node )) {
      print join("\t", $node, $edge, $graph{$node}[ $OUT ]{$edge}, join("-", sort {$a <=> $b} keys %{_reads($node)}),"\n");
#      print join("\t", $node, $edge, $graph{$node}[ $OUT ]{$edge}, int(keys %{_reads($node)}),int(keys %{_reads($edge)}),"\n");
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
      if ( $seq && length($seq) > 10) {
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

#  $REVERSE = 1 if ($filename =~ /R/);
  
  my ($reads, $kmer_counts);# = count_kmers( $filename );

  my ($name, $seq) = (1, "");
  my $exit_counter = $EXIT_COUNTER;

  $infile = $filename;

  open( my $in, $filename) || die "Could not open '$filename': $!\n";
  while(<$in>) {
    chomp;

    if ( /\>/  ) {
      if ( $seq && length($seq) >= $FRAG_MIN && length($seq) <= $FRAG_MAX ) {
#      print "$name\n";
	DeBruijn::add_sequence( $name, $seq, $kmer_counts, $reads );
      }
      $seq = "";
      last if ( $exit_counter-- == 0 );
      $name++;
    }
    else {
      $seq .= $_;
    }
  }

  DeBruijn::add_sequence( $name, $seq, $kmer_counts, $reads ) if ( $seq && length($seq) >= $FRAG_MIN && length($seq) <= $FRAG_MAX);

  # uses the multipliers table to adjust the weights of the nodes.
#  _adjust_weights();
#x  warn qx{ ps -o rss,vsz $$ }, "\n";
#  print Dumper( \%graph);
}



# 
# 
# 
# Kim Brugger (28 Oct 2011)
sub _adjust_weights {
  
  dump_graph();

#  foreach my $node ( keys %graph ) {
#    my $reads = keys %{_reads( $node )};


}



# 
# 
# 
# Kim Brugger (19 Oct 2011)
sub cyclic_sequence {
  my ( $seq ) = @_;

  my %kmers;
  for( my $i = 0; $i<length($seq) - $kmer; $i++) {
    my $node = substr( $seq, $i, $kmer);
    
    $kmers{ $node }++;
#    print "repeated sequence: $node \n" if ( $kmers{ $node } > 1);
    
    if ( $kmers{ $node } > 1) {
#      $seq =~ s/(.*?$node.*?)$node.*/$1/;
#      print "Cleaned sequence: $seq\n";
      return 1;
    }
  }
  return 0;
}


# 
# 
# 
# Kim Brugger (06 Oct 2011)
sub add_sequence {
  my ($name, $seq, $kmer_counts, $reads) = @_;

  $seq = reverse($seq) if ($REVERSE);
  
#  print "$name\t$seq\n";

  if ( cyclic_sequence( $seq ) ) {
#    print  "Ignoring cyclic sequence ($name, $seq)\n";
    return;
  }

#  print "$seq\n";

  $build_reads++;
  my $prev_node = "S";

  if (0 && $multipliers{ $seq } ) {
    $multipliers{ $seq }++;
    return;
  }

  $multipliers{ $seq }++;

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

  my $start = 'S';

#  print "start paths: ". join(" ", keys %{$graph{'S'}{'OUT'}}) . "\n";

#  print Dumper(\%graph);

  foreach my $start_pos ( _out( 'S' )) {
    my @post_poss = _out( $start_pos );

    if ( ! @post_poss ) {
      my $weight = 1;#keys %$legacy;
      $start_pos = reverse($start_pos) if ($REVERSE);

      print ">$infile\n$start_pos\n" if ( length($start_pos) >= $FRAG_MIN );
      next;
    }

    my $start_reads = _reads( $start_pos );
#    print Dumper( $start_reads );
    foreach my $post_pos ( @post_poss ) {

#      print "S $start_pos $post_pos\n";
      my $post_reads = _reads( $post_pos );
#      print Dumper( $post_reads );

      my %shared_reads;
      foreach my $read ( keys %$post_reads ) {
	if ( $$start_reads{ $read }) {
#	  print "SHARED READ S:: $start_pos -- $post_pos $read :: \n";
	  $shared_reads{ $read } = ($$start_reads{ $read } || 0) + ($$post_reads{ $read } || 0);
	}
      }

      if ( keys %shared_reads ) {
	_path_finder($start_pos.substr($post_pos, $kmer -1), $post_pos, \%shared_reads);
      }
    }
  }
}



# 
# 
# 
# Kim Brugger (07 Oct 2011)
sub _path_finder {
  my ($pre_path, $pos, $prev_shared_reads) = @_;

#  print "PREV " . Dumper( $prev_shared_reads );

#  print "$pre_path -- $pos\n";


  my @post_poss = _out( $pos );
  if (! @post_poss ) {
    my $weight = 1;#keys %$legacy;
    $pre_path = reverse($pre_path) if ($REVERSE);

    print ">$infile\n$pre_path\n" if ( length($pre_path) >= $FRAG_MIN );
    return;
  }



  my %shared_reads;

  foreach my $post_pos ( @post_poss ) {
#    print "Trying to connect $pos with $post_pos\n";

    my $post_reads = _reads( $post_pos );
#    print Dumper( $post_reads );

    foreach my $read ( keys %$post_reads ) {
      if ( $$prev_shared_reads{ $read }) {
#	print "SHARED READ S:: $pos -- $post_pos $read :: \n";
	$shared_reads{ $post_pos}{$read } = $$prev_shared_reads{$read} + $$post_reads{$read};
      }
    }
  }
  
 
#  print Dumper( \%shared_reads );
 
  if ( keys %shared_reads > 1) {
    
    my %ratings;
    foreach my $outgoing ( keys %shared_reads ) {
      foreach my $read ( keys %{$shared_reads{ $outgoing }} ) {
	$ratings{$outgoing} += $shared_reads{ $outgoing }{ $read };
      }
    }
    
    print Dumper(\%ratings);

    my $outgoing = (sort { $ratings{$b} <=> $ratings{$a}} keys %ratings)[0];
    print "Picked $outgoing\n";
    my $out_reads = _reads( $outgoing );
    
    foreach my $read ( keys %$out_reads ) {
#      $$prev_shared_reads{ $read } ||= 0;
#      $$prev_shared_reads{ $read } += $$out_reads{$read};
      $$prev_shared_reads{ $read } += $$out_reads{ $read } if ( $$prev_shared_reads{ $read });
    }
    _path_finder ($pre_path .substr($outgoing, $kmer -1), $outgoing, $prev_shared_reads );

    ;
  }
  elsif ( keys %shared_reads == 1) {
#    return;
#    print "SINGLE SHARED READ, WTF!: ".( join(", ", keys %shared_reads ))."\n";
    my $outgoing = (keys %shared_reads)[0];
    my $out_reads = _reads( $outgoing );

    foreach my $read ( keys %$out_reads ) {
      $$prev_shared_reads{ $read } ||= 0;
      $$prev_shared_reads{ $read } += $$out_reads{ $read } if ( $$prev_shared_reads{ $read });
    }
    _path_finder ($pre_path .substr($outgoing, $kmer -1), $outgoing, $prev_shared_reads );
  }
  else {
    print Dumper( \%shared_reads );
    print "No shared reads from $pos onwards\n";
  }
  
}


1;



__END__



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


