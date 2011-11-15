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

#$file =~ s/.fa/.tkr/;
#Typer::dump_db($file)  if ($file =~ /.tkr/);
Typer::readin_db($file)  if ($file =~ /.tkr/);

my($name, $seq);
while(<>) {
  chomp;
    
  if ( /\>/  ) {
    if ( $seq ) {
      print "$name\n";
      my ($score, $types) = Typer::type_sequence( $seq );
#      @{$types} = grep(/0[37]:01/, @{$types});
      my $type_head = substr(join(" ", @{$types}), 0 , 50);
      $type_head = join(" ", sort @{$types});
      
      print "SCORE :: $score, $type_head\n";
    }
    $seq = "";
    $name = $_;
  }
  else {
    $seq .= $_;
  }
}
if ( $name && $seq ) {
      print "$name\n";
      my ($score, $types) = Typer::type_sequence( $seq );
#      @{$types} = grep(/0[37]:01/, @{$types});
      my $type_head = substr(join(" ", @{$types}), 0 , 50);
      $type_head = join(" ", sort @{$types});
      print "SCORE :: $score, $type_head\n";
}
  
