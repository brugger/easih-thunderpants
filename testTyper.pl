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

Typer::readin_file(shift);
Typer::type_sequence("CACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTACCTGGACAGATACTTCCATAACCAGGAGGAGAACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGGCCGGGTGGACAACTACTGCAGACACAACTACGGGGTTGTGGAGAGCTTCACAGTGCAGCGGCGAG");

Typer::type_sequence("GCGGCGTCGCTGTCAGTGTCTTCTCAGGAGGCCGCCCGTGTGACCGGATCCTTCGTGTACCCGCAGCACGTTTCTTGGAGCTGCGTAAGTCTGAGTGTCATTT");
Typer::type_sequence("GCGGCGTCACTGTCAGTGTCTTCTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTT");
Typer::type_sequence("GCGCGTCGCTGTCAGTGTTTTTCCCGGAGACCGCCCCTGTGACCGGATCGTTTGTGCCCCCACAGCACGTTTCCTGTGGCAGGGTAAGTATAAGTGTCATTT");
Typer::type_sequence("GCGCGTCGCTGTCAGTGTTTTTCCCGGAGACCGCCCCTGTGACCGGATCGTTTGTGCCCCACAGCACGTTTCCTGTGGCAGGGTAAGTATAAGTGTCATTT");
Typer::type_sequence("GCGGCGTCGCTGTCAGTGTTTTTCCCGGAGACCGCCCCTGTGACCGGATCGTTTGTGCCCCCACAGCACGTTTCCTGTGGCAGGGTAAGTATAAGTGTCATTT");
Typer::type_sequence("GCGGCGTCGCTGTCAGTGTTTTCCCGGAGACCGCCCCTGTGACCGGATCGTTTGTGCCCCCACAGCACGTTTCCTGTGGCAGGGTAAGTATAAGTGTCATTT");
Typer::type_sequence("GCGGCGTCGCTGTCAGTGTCTTCTCAGGAGGCCGCCCGTGTGACCGGATCCTTCGTGTACCCGCAGCACGTTTCCTGTGGCAGGGTAAGTATAAGTGTCATTT");
Typer::type_sequence("GCGGCGTCACTGTCAGTGTCTTCTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTT");
Typer::type_sequence("GCGGCGTCGCTGTCAGTGTCTTCTCAGGAGGCCGCCCGTGTGACCGGATCCTTCGTGTACCCGCAGCACGTTTCTTGGAGCTGCGTAAGTCTGAGTGTCATTT");
Typer::type_sequence("CACGTTTCCTGTGGCAGGGTAAGTATAAGTGTCATTTCTTCAACGGGACGGAGCGGGTGCAGTTCCTGGAAAGACTCTTCTATAACCAGGAGGAGTTCGTGCGCTTCG");
Typer::type_sequence("CACGTTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTACCTGGACAGATACTTCCATAACCAGGAGGAGAACGTGCGCTTCG");
Typer::type_sequence("CACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTACCTGGACAGATACTTCCATAACCAGGAGGAGAACGTGCGCTTCG");
Typer::type_sequence("CACGTTTCCTGTGGCAGGGTAAGTATAAGTGTCATTTCTTCAACGGGACGGAGCGGGTGCAGTTCCTGGAAAGACTCTTCTATAACCAGGAGGAGTTCGTGCGCTTCG");
Typer::type_sequence("CACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTACCTGGACAGATACTTCCATAACCAGGAGGAGAACGTGCGCTTCG");
Typer::type_sequence("CGGGTGCAGTTCCTGGAAAGACTCTTCTATAACCAGGAGGAGTTCGTGCGCTTCGACAGCGACGTGGGGGAGTACCGGGCGGTGACGGAGCTAGGGCGGCCTGTCGCCGAGTCCTGG");
Typer::type_sequence("CGGGTGCGGTACCTGGACAGATACTTCCATAACCAGGAGGAGAACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTG");
Typer::type_sequence("CGGGTGCGGTACCTGGACAGATACTTCCATAACCAGGAGGAGTTCCTGCGCTTCGACAGCGACGTGGGGGAGTACCGGGCGGTGACGGAGCTGGGGCGGCCTGTCGCCGAGTCCTG");
Typer::type_sequence("GGGGGAGTACCGGGCGGTGACGGAGCTAGGGCGGCCTGTCGCCGAGTCCTGGAACAGCCAGAAGGACATCCTGGAGGACAGGCGGGGCCAGGTGGACACCGTG");
Typer::type_sequence("GGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGGCCGGGTGGACAACTAC");
Typer::type_sequence("GGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGGCCGGGTGGACAACTAC");
Typer::type_sequence("GGGGGAGTACCGGGCGGTGACGGAGCTAGGGCGGCCTGTCGCCGAGTCCTGGAACAGCCAGAAGGACATCCTGGAGGACAGGCGGGCCAGGTGGACACCGTG");
Typer::type_sequence("CTCCTGGAGCAGAAGCGGGGCCGGGTGGACAACTACTGCAGACACAACTACGGGGTTGTGGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGCGGGGCGGGGC");
Typer::type_sequence("CTCCTGGAGCAGAAGCGGGGCCGGGTGGACAACTACTGCAGACACAACTACGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCATGTCGGGGGGCGGGGC");
Typer::type_sequence("CTCCTGGAGCAGAAGCGGGGCCGGGTGGACAATTACTGCAGACACAACTACGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCATGTCGGGGGGCGGGGC");
Typer::type_sequence("ATCCTGGAGGACAGGCGGGGCCAGGTGGACACCGTGTGCAGACACAACTACGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGGAGGCGGGGC");
Typer::type_sequence("ATCCTGGAGGACAGGCGGGGCCAGGTGGACACCGTGTGCAGACACAACTACGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCATGTCGGGGGGCGGGGC");
Typer::type_sequence("CTCCTGGAGCAGAAGCGGGGCCGGGTGGACAACTACTGCAGACACAACTACGGGGTTGTGGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGCGGGGCGGGGC");
Typer::type_sequence("CTCCTGGAGCAGAAGCGGGCCGGGTGGACAATTACTGCAGACACAACTACGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCATGTCGGGGGCGGGC");
Typer::type_sequence("ATCCTGGAGGACAGGCGGGCCAGGTGGACACCGTGTGCAGACACAACTACGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGGAGGCGGGGC");

