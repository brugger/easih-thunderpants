#!/usr/bin/perl -w
# 
# Recalibrates QV values by calculating the true error rate for all QVs.
# 
#
# Kim Brugger (11 Oct 2011), contact: kim.brugger@easih.ac.uk

use strict;
use Data::Dumper;
use warnings;
use Getopt::Std;

my $PRE_POS     = 1;
my $PRE_PRE_POS = 0;
my $MIN_SPLIT      = $opts{'s'} || 60; # should probably be something like 10-15...
my $MIN_DEPTH      = $opts{'d'} || 15; # I guess something between 7 and 20 should be good


my %opts;
getopts('b:d:R:B:hs:Sr:UM:m:L:l:s:', \%opts);
my $bam_file = $opts{'b'} || Usage();
my $chr_file = $opts{'R'} || Usage();

die "'$chr_file' does not exist\n" 
    if ( ! -e $chr_file );
die "index does not exist for '$chr_file', please create one with samtools faidx\n"
    if (! -e "$chr_file.fai");

my $samtools    = find_program('samtools');

my $BWA = 1;


my $FILTER_BUFFER  = $opts{'B'} || 100;
my $sam_out        = $opts{'S'};

my $log_file       = $opts{'L'} || 0;
my $report_file    = $opts{'l'};


open (my $log, "> $log_file") || die "Could not open file '$log_file': $!\n" if ( $log_file );

my $good_trans = legal_transitions();

open (STDOUT, " | $samtools view -Sb - ") || die "Could not open bam stream: $!\n" if ( !$sam_out );
print STDOUT  `$samtools view -H $bam_file`;
my ( $s, $d, $e, $b_pre, $b_post) = (0,0,0, 0, 0);

my $region      = $opts{'r'};
my $sample_size = $opts{'s'} || 0;
$region = "chr:200-400";

my $fasta;

my %QV_stats;

if ( ! $region ) {
  # loop through the regions one by one, but only keep the most current chr in memory.
  
  open(my $spipe, "$samtools view -H $bam_file | ") || die "Could not open '$bam_file': $!\n";
  while(<$spipe>) {
    next if ( ! /\@SQ/);
    foreach my $field ( split("\t") ) {
      if ( $field =~ /SN:(.*)/) {
	$fasta = readfasta( $chr_file, $1 );
	analyse($1)
      }
    }
  }
}
else {
  $region =~ s/,//g;
  if ($region =~ /^(\w+):\d+-\d+/ || $region =~ /^(\w+):\d+\z/ ) {
    $fasta = readfasta( $chr_file, $1 );
    analyse($region)
  }
  else {
    $fasta = readfasta( $chr_file, $region );
    analyse($region)
  }
}


if ($report_file) {
  open( my $report, "> $report_file") || die "Could not write to '$report_file': $!\n";
  print $report "Scrubbed ". ( $s+$d+$e) . " reads from $bam_file\n";
  print $report "\t$s single colour errors\n";
  print $report "\t$d double colour errors\n";
  print $report "\t$e end    colour errors\n";
}

# 
# 
# 
# Kim Brugger (05 Nov 2010)
sub analyse {
  my ($region) = @_;

  my (@splits, @ref_ids, @reads, @SNPs);

  my $current_pos = undef;

  open (my $bam, "$samtools view $bam_file $region | ") || die "Could not open 'stream': $!\n";
  while(<$bam>) {

    chomp;
    my @F = split("\t");
    my ($id, $flags, $chr, $pos, $mapq, $cigar, $mate, $mate_pos, $insert_size, $seq, $quality, @opts) = @F;
    
    my $entry = { sam    => \@F,
		  id     => $id,
		  flags  => $flags, 
		  chr    => $chr, 
		  cigar  => $cigar,
		  pos    => $pos};
    

    # Not a mapped read
    if ($flags & 0x0004 ) {
      next;
    }

    if ( 0 && $current_pos && $current_pos != $pos )  {
      # the cs_splits array needs to be synced with this new pos. Bring forth the 
      # array as many places. If there is a gap, traverse the whole thing and reset
      # the whole thing, so we start from fresh.
      
      print STDERR "stepping from $current_pos =>>  $pos\n";
      for(my $i = 0; $i < $pos - $current_pos; $i++ ) {

	# This is a sliding array, keeping track of all the colour balances, and update the SNP array...
	shift @ref_ids if ( @ref_ids >= 2 );
	
	next if (! $splits[0]{total} || $splits[0]{total} == 0);
	
	
	my ($right) = (0);
	map { $right += $splits[0]{ $_ } if ( $splits[0]{ $_ } && $splits[0]{ref} eq $_)} ( 'O','1','2','3');
	
	push @ref_ids, [$splits[0]{pos}, $right*100/$splits[0]{total}, $splits[0]{total}];
	shift @splits;
	
	if ( @ref_ids == 2 && 
	     $ref_ids[ $PRE_PRE_POS ] && $ref_ids[ $PRE_POS ]  &&
	     $ref_ids[ $PRE_PRE_POS ][ 1 ] < $MIN_SPLIT && 
	     $ref_ids[ $PRE_POS ][ 1 ]     < $MIN_SPLIT && 
	     $ref_ids[ $PRE_POS ][ 2 ]     > $MIN_DEPTH     &&
	     $ref_ids[ $PRE_PRE_POS ][ 2 ] > $MIN_DEPTH ) {


	  print STDERR "SNP at pos $region:".($current_pos  + $i)." (depth: $ref_ids[ $PRE_POS ][ 2 ]) splits: ($ref_ids[ $PRE_PRE_POS ][1], $ref_ids[ $PRE_POS ][1])\n";
	  push @SNPs, [$ref_ids[ $PRE_PRE_POS ], $ref_ids[ $PRE_POS ]];
	}
	
	if ( @reads &&  $reads[0]{ pos } + $FILTER_BUFFER <= $pos + $i ) {
	  #remove SNPs further behind that we will ever be seeing again...
#	print STDERR "Buffered SNPs: " . @SNPs . "\n". " ($reads[0]{ pos } + $FILTER_BUFFER > $SNPs[0][1][0]); \n";
	  while ( @SNPs ) { 
#	  print STDERR Dumper( $SNPs[0]);
	    
	    last 	if ($reads[0]{ pos } - $FILTER_BUFFER < $SNPs[0][1][0]);
	    my @snp = shift @SNPs;
	  }
#	print STDERR "Buffered SNPs: " . @SNPs . "\n". " ($reads[0]{ pos } + $FILTER_BUFFER > $SNPs[0][1][0]); \n";
#	print STDERR "---------------------------------\n";
#	  exit;
	  
	  #scrub all the reads that are far enough away..
	  while ( @reads > 0 && $reads[0]{ pos } + $FILTER_BUFFER < $pos + $i ) {
	    my $read = shift @reads;
	    
	    scrub( $read, \@SNPs);
	    print_sam( $read );
	  }
	  
#	print "Post-removal nr of reads: " . @reads ."\n";
	  
	}
	
      }
    }
    $current_pos = $pos;
    
    
#  print "pos: $current_pos, buffered reads: " . @reads . "\n";
    
    
    my $hlen = length($csfasta);
    if ( 1 &&  $cigar =~ /[IDNP]/) {
      my $t_cigar = $cigar;
      my $padding = 2;
      $t_cigar =~ s/(\d+)[ID]/ {$padding += $1}/ge;
      $hlen += $padding;
    }

    my $gseq = substr($fasta, $pos - 1, $hlen);

    # for some odd reason, I cannot be arsed to figure out right now...
#  $gseq = substr($gseq, 0, length($csfasta)) if ( length( $gseq) > length($csfasta));
    
    # pseq: patched sequence, the aligned sequence
    # pqual: the quality values for the aligned bases.
    my ($pseq, $pqual) = patch_alignment($seq, $qual, $cigar);

    $$entry{end  } = $pos + length($pseq);
    $$entry{pseq } = $pseq;
    $$entry{pqual} = $pqual;

    push @reads, $entry;
#  print "$csfasta\n";
    
    my @gseq = split("", $gseq);
    my @pseq = split("", $pseq);
    my $gaps = 0;
    for(my $i = 0; $i<@gesq; $i++ ) {
      
      # there is an insert in the reference, so this number needs to 
      # be subtracted to get the real genome position.
      if ($gcsf[ $i ] eq "-") {
	$gaps++;
	next;
      }

      $splits[$i - $gaps]{ ref      } = $gseq[ $i ];
      $splits[$i - $gaps]{ pos      } = $pos  + $i;
      $splits[$i - $gaps]{ B  }{ $pseq[$i]  }++;
      $splits[$i - $gaps]{ QV }{ $pqual[$i] }++;
      next if ($pseq[ $i ] eq "-");
      $splits[$i - $gaps]{ total    }++;
      
    }
  }


  print Dumper( \%splits );

# empty the reads buffer..
  while ( @reads ) {
    my $read = shift @reads;
    
    scrub( $read, \@SNPs);
    print_sam( $read );
    
  }
}




# 
# 
# 
# Kim Brugger (28 Oct 2010)
sub print_sam {
  my ( $read ) = @_;

#  return;
  my $sam    = $$read{ sam   };
  @$sam[ 1 ] = $$read{ flags };
  @$sam[ 4 ] = $$read{ mapq  };
  
  print STDOUT join("\t", @$sam) . "\n";
}


# 
# 
# 
# Kim Brugger (28 Oct 2010)
sub scrub {
  my ( $read, $SNPs ) = @_;

  return if ( $$read{indel} );

  return if ( ! $$read{singles} && ! $$read{doubles} );
	  
  foreach my $snp ( @$SNPs ) {
    my ($snp_start, $snp_end) = ($$snp[0][0], $$snp[1][0]);
#    next if ( ! $snp_start || ! $snp_end);

    if ( $$read{ pos } - 5 > $snp_end) {
#      print STDERR "Pre-Bailing... $$read{pos} -> $$read{end} vs $snp_start => $snp_end\n";
      $b_post++;
      next;
    }


    if ( $$read{ end } + 5 < $snp_start) {
#      print STDERR "Bailing...\n";
      $b_pre++;
      return;
    }

#    next;
      
    foreach my $single (@{$$read{singles}}) {

      if (($single - 1 + $$read{pos} <= $snp_start &&
	   $single + 1 + $$read{pos} >= $snp_end)
	  || 
	  ($single - 1 + $$read{pos} <= $snp_start &&
	   $single + 1 + $$read{pos} >= $snp_start)
	  ||
	  ($single - 1 + $$read{pos} <= $snp_end &&
	   $single + 1 + $$read{pos} >= $snp_end)) {

	print STDERR "$$read{id} -- $$read{cigar} $$read{flags}\n$$read{a}";
	
	
	$$read{flags} += 4 if ( $set_unmapped );
	$$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);
	push @{$$read{sam}}, "ES:Z:1";
	$s++;
	return;
      }
    }  

#    next;

    foreach my $double (@{$$read{doubles}}) {
      
      if (($$double[0] + $$read{pos} <= $snp_start &&
	   $$double[1] + $$read{pos} >= $snp_end)
	  || 
	  ($$double[0] + $$read{pos} <= $snp_start &&
	   $$double[1] + $$read{pos} >= $snp_start)
	  ||
	  ($$double[0] + $$read{pos} <= $snp_end &&
	   $$double[1] + $$read{pos} >= $snp_end)) {
	
#	print STDERR "$$read{id} -- $$read{cigar}\n$$read{a}";
#	print $$read{a};

	$$read{flags} += 4 if ( $set_unmapped );
	$$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);
	push @{$$read{sam}}, "ES:Z:2";
	$d++;
	return;
      }
    }  	  	  
	  
    if ( $$read{pos} == $snp_end || $$read{end}  == $snp_end ) {
	    
#      print STDERR "$$read{id} -- $$read{cigar}\n$$read{a}";
#      print $$read{a};
      $$read{flags} -= 4 if ( $set_unmapped );
      $$read{mapq}   = $set_mapq_score if ( defined $set_mapq_score);

      push @{$$read{sam}}, "ES:Z:3";
#      $$read{flags} -= 4;
#      $$read{mapq}   = 2;		
      $e++;
      return;
    }
  }

  return;
}



# 
# probe vs genome
# 
# Kim Brugger (13 Oct 2010)
sub align {
  my ( $s1, $s2, $strand) = @_;

  return ([],[], "NA\n") if ($s1 eq $s2);

  my @s1 = split("", $s1);
  my @s2 = split("", $s2);

  my (@singles, @align, @doubles);
  my $snp = 0;
  
  for(my $i = 0; $i < @s1; $i++) {
    if ( $s1[$i] eq $s2[$i] ) {
      $align[ $i ] = 1;
      $snp = 0;
    }
    else {
      $align[ $i ] = 0;
      # Check and see if this is a double "error" IE a SNP
      if ( $i > 0 && ! $align[ $i - 1 ]) {

	if ( !$snp ) {
	  push @doubles, [$singles[-1], $i];
	  pop @singles;

	}
	$snp = 1;
	$doubles[-1][1] = $i;
      }
      else {
	push @singles, $i;
	$snp = 0;
      }
    }
  }

  my @long;

  foreach my $d ( @doubles ) {
    
    my $s1_d = join("",@s1[$$d[0]..$$d[1]]);
    my $s2_d = join("",@s2[$$d[0]..$$d[1]]);

    if ( ! $$good_trans{$s1_d}{$s2_d}) {
      push @long, $d;      
    }
    
  }

  if ( 1 ) {
    my $align = join("", @align);
    $align =~ tr/01/ \|/;
#    print STDERR "$s1\n$align\n$s2\n\n" if ( @singles || @long);
    return( \@singles, \@long, "$s1\n$align\n$s2\n\n");
  }

  return( \@singles, \@long, );
}






sub patch_alignment {
  my ( $read, $ref, $cigar ) = @_;

  # Extended cigar format definition ( from the sam/bam format file)
  # M Alignment match (can be a sequence match or mismatch)
  # I Insertion to the reference
  # D Deletion from the reference
  # N Skipped region from the reference
  # S Soft clip on the read (clipped sequence present in <seq>)
  # H Hard clip on the read (clipped sequence NOT present in <seq>)
  # P Padding (silent deletion from the padded reference sequence)

  if ( $cigar !~ /[HDIS]/) {
    $ref = substr($ref, 0, length($read));
    return ($read, $ref);
  }
  
  my @read  = split("", $read );
  my @ref   = split("", $ref );
  
  my $ref_cigar = $cigar;
  $ref_cigar =~ s/^\d+[HDS]//;

  my (@cigar) = $ref_cigar =~ /(\d+\w)/g;

  my $offset = 0;
  foreach my $patch ( @cigar ) {
    my ($length, $type) =  $patch =~ /(\d+)(\w)/;

    if ( $type eq 'M') {
      $offset += $length;
      next;
    }
    elsif ( $type eq "I") {
      my @dashes = split("", "-"x$length);
      splice(@ref, $offset, 0, @dashes);
    }
    elsif ( $type eq "S" || $type eq "H") {
      splice(@ref,  $offset, $length);
    }    
  }

  if (1){
    $offset = 0;
    my (@cigar) = $cigar =~ /(\d+\w)/g;
    foreach my $patch ( @cigar ) {
      my ($length, $type) =  $patch =~ /^(\d+)(\w)/;

      if ( $type eq 'M') {
	$offset += $length;
	next;
      }
      elsif ( $type eq "D") {
	my @dashes = split("", "-"x$length);
	splice(@read,  $offset, 0, @dashes);
	$offset += $length;
      }
      elsif ( $type eq "H" || $type eq "S") {
      splice(@read,  $offset, $length);
      }    
    
    }
  }
  $read = join("", @read);
  $ref  = join("", @ref );

  $ref = substr($ref, 0, length($read));
  
  return ($read, $ref);
}




#
# Read the fasta files and puts entries into a nice array
#
sub readfasta {
  my ($file, $region) = @_;  

  my $sequence;
  my $header;
  
  open (my $f, "$samtools faidx $file $region |" ) || die "Could not open $file:$1\n";
  while (<$f>) {
    chomp;
    if (/^\>/) {
      if ($header) { # we have a name and a seq
	return ($header, $sequence);
      }
      $header = $_;
      $header =~ s/^\>//;
    }
    else {$sequence .= $_;}
  }
  

  return ($sequence);
}




# 
# 
# 
# Kim Brugger (13 Jul 2010)
sub find_program {
  my ($program) = @_;
  
  my $username = scalar getpwuid $<;
  
  my @paths = ("/home/$username/bin/",
	       "./",
	       "/usr/local/bin");
  
  foreach my $path ( @paths ) {
    return "$path/$program" if ( -e "$path/$program" );
  }

  my $location = `which $program`;
  chomp( $location);
  
  return $location if ( $location );
  
  die "Could not find '$program'\n";
}


# 
# 
# 
# Kim Brugger (05 Nov 2010)
sub Usage {
  $0 =~ s/.*\///;
  die "USAGE: $0 -b<am file> -R<eference genome (fasta)> -d[ min depth, default=15] -s[ min Split, default=60] -B[uffer, default=100] -M[ set mapq score for offending reads] -U[n set mapped flag for offending reads]\n";

  # tests :::: odd SNP reporting:  10:74879852-74879852
  # large indel: 10:111800742
}
