package Typer;
# 
# 
# 
# 
# Kim Brugger (04 Nov 2011), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
no warnings "recursion";
use Storable;

my $kmer = 12;
my %kmers;

my $EXIT_COUNTER = -1;


# 
# 
# 
# Kim Brugger (06 Oct 2011)
sub type_sequence {
  my ($seq) = @_;
  
  my (%dna_hits, %aa_hits);

  my ($steps, $aa_steps) = (0,0);

  for( my $i = 0; $i<length($seq) - $kmer+1; $i++) {
    my $node = substr( $seq, $i, $kmer);
    my $AAnode = translate($node);
    
    map { $dna_hits{ $_ }++ } keys %{$kmers{ $node }};
    $steps++;

    if ( $AAnode !~ /\./) {
      map {  $aa_hits{ $_ }++ } keys %{$kmers{ $AAnode }};
      $aa_steps++;
    }
  }

  #pull out the best hit,DNA first.
  my (%weighted_dna_types, $dna_heaviest);
  map { push @{$weighted_dna_types{$dna_hits{$_}}}, $_;
	$dna_heaviest = $dna_hits{$_} if ( !$dna_heaviest || $dna_heaviest < $dna_hits{$_});
  } sort keys %dna_hits;

  my (%weighted_aa_types, $aa_heaviest);
  map { push @{$weighted_aa_types{$aa_hits{$_}}}, $_;
	$aa_heaviest = $aa_hits{$_} if ( !$aa_heaviest || $aa_heaviest < $aa_hits{$_});
  } sort keys %aa_hits;

#  my $dna_string = join(" ", @{$weighted_dna_types{ $dna_heaviest }});
#  my $aa_string = join(" ", @{$weighted_aa_types{ $aa_heaviest }});

#  return if ($dna_heaviest != $steps && $aa_heaviest != $aa_steps);
  print "Seq splits into $steps/$aa_steps DNA/AA fragments. ($dna_heaviest/$aa_heaviest)\n";

#  print "=== DNA & AA agrees ===\n" if ( $dna_string eq $aa_string);

#  print substr($dna_string, 0, 20) . "\n" if ( $dna_string eq $aa_string);

  return (100, $weighted_dna_types{ $dna_heaviest }) if ( $dna_heaviest == $steps);
  return (99 , $weighted_aa_types{ $aa_heaviest })   if ( $aa_heaviest == $aa_steps);
#  return ( 99, $weighted_dna_types{ $dna_heaviest }) if ( $dna_string eq $aa_string);

  if ( $steps - $dna_heaviest > $aa_steps - $aa_heaviest) {
    return (sprintf("%.2f", 100*$dna_heaviest/$steps), $weighted_aa_types{ $aa_heaviest });
  }
  else {
#    return (sprintf("%.2f", 100*($aa_heaviest/$aa_steps)), $weighted_aa_types{ $aa_heaviest });
  }

  return (0, []);
}


# 
# Translate a DNA string in all 6 frames
# 
# Kim Brugger (05 Feb 2004)
sub six_frames {
  my ($DNA) = @_;

  my $res = "";

  for (my $i = 0; $i < 3; $i++) {
    $res .= ">+strand frame $i\n";

    my $DNAseq = substr ($DNA, $i);
    $res .= nicefasta(translate($DNAseq, 1));
  }

  my $revDNA = revDNA($DNA);
  for (my $i = 0; $i < 3; $i++) {
    $res .= ">-strand frame $i\n";

    my $revDNAseq = substr ($revDNA, $i);
    $res .= nicefasta(translate($revDNAseq, 1));
  }

  return $res;
}


sub translate {
  my($sequence, 
     $multi_M, # Uses M where possible (RAW translate ?)
     ) = (@_);

  $sequence =~ tr/[atgc]/[ATGC]/;
  my %alt_start_codon = ("ATG" => "M", "GTG" => "M", "TTG" => "M");




  my %code = ("ATG" => "M", 
              "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
              "CGT" => "R", "CGC" => "R", "CGA" => "R", 
              "CGG" => "R", "AGA" => "R", "AGG" => "R",
              "AAT" => "N", "AAC" => "N", 
              "GAT" => "D", "GAC" => "D", "TGT" => "C",
              "TGC" => "C",
              "CAA" => "Q", "CAG" => "Q",
              "GAA" => "E", "GAG" => "E",
              "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", 
              "CAT" => "H", "CAC" => "H",
              "ATT" => "I", "ATC" => "I", "ATA" => "I",
              "TTA" => "L", "CTT" => "L", "TTG" => "L",
              "CTC" => "L", "CTA" => "L", "CTG" => "L", 
              "AAG" => "K", "AAA" => "K", 
              "TTT" => "F", "TTC" => "F",
              "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P", 
              "AGT" => "S", "AGC" => "S", "TCT" => "S", 
              "TCC" => "S", "TCA" => "S", "TCG" => "S", 
              "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T", 
              "TGG" => "W",
              "TAT" => "Y", "TAC" => "Y",
              "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V", 
              "TAG" => ".", "TAA" => ".", "TGA" => ".");
  
  my $translated = "";
  my $first_M = 1;

  # The last 3 N gets translated into a start codon so we can be sure
  # that partial genes are also identified
#  $sequence =~ s/(N{20,})NNN/$1ATG/g;
  for (my $i = 0; $i < length $sequence; $i+=3) {
    my $codon = substr($sequence,$i,3);
    $codon =~ tr/[atgc]/[ATGC]/;
    if ($code{$codon}) {
      if ( $alt_start_codon{$codon} && $first_M ) {
        $translated .= $alt_start_codon{$codon};
        $first_M = 0;
      }
      elsif ($multi_M && $alt_start_codon{$codon}) {
        $translated .= $alt_start_codon{$codon};
      }
      else {
        $translated .= $code{$codon};
      }
    }
    elsif ($codon =~ /X/) {
      $translated .= ".";
    }
    elsif ($codon =~ /NNN/) {
      $translated .= "X";
    }
    # The last N(S) gets translated into a start codon so we can be sure
    # that partial genes are also identified
    elsif ($codon =~ /N[ATGC]{1,2}/) {
      $translated .= "M";
    }
    # We do not know how to translate this so insert an X
    else {
      $translated .= "X";
    }
  }

  return $translated;
}


# 
# 
# 
# Kim Brugger (14 Oct 2011)
sub set_kmer {
  my ($new_kmer) = @_;

  if ($new_kmer !~ /^\d+\z/ || 
      $new_kmer <= 0 ||
      $new_kmer >  53 ) {
    
    print STDERR "kmers should be a integer between 0 and 53 not '$new_kmer'\n";
    return;
  }

  $kmer = $new_kmer;
}






# 
# 
# 
# Kim Brugger (09 Nov 2011)
sub dump_db {
  my ($outfile) = @_;

  die "Need an outfile to write the database to\n" if (! $outfile );

  return Storable::store(\%kmers, $outfile);
}



# 
# 
# 
# Kim Brugger (09 Nov 2011)
sub readin_db {
  my ($dbfile) = @_;

  die "Need an dbfile to read from\n" if (!$dbfile);


  %kmers = %{Storable::retrieve( $dbfile)};
}



# 
# 
# 
# Kim Brugger (17 Oct 2011)
sub readin_file {
  my ($filename) = @_;

  
  my ($reads, $kmer_counts);# = count_kmers( $filename );

  my ($name, $seq);
  my $exit_counter = $EXIT_COUNTER;



  open( my $in, $filename) || die "Could not open '$filename': $!\n";
  while(<$in>) {
    chomp;

    if ( /\>/  ) {
      
      _add_sequence( $name, $seq ) if ( $seq );
      $seq = "";
      $name = $_;
      $name =~ s/\>//;
      last if ( $exit_counter-- == 0 );
    }
    else {
      $seq .= $_;
    }
  }

  _add_sequence( $name, $seq ) if ( $seq );
#  print Dumper( \%kmers );
#  exit;
}


# 
# 
# 
# Kim Brugger (06 Oct 2011)
sub _add_sequence {
  my ($name, $seq) = @_;

#  print "$name\n";
#  print "$seq\n";

  for( my $i = 0; $i<length($seq) - $kmer+1; $i++) {
    my $node = substr( $seq, $i, $kmer);

    my $AAnode = translate($node);
    $kmers{ $node }{ $name }++;
    $kmers{ $AAnode }{ $name }++ if ( $AAnode !~ /\./);
  }

}


1;



__END__


