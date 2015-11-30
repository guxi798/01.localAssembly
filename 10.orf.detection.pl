#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $seqfile = shift @ARGV;
my $format = shift @ARGV;

my $seqio_obj = Bio::SeqIO->new(-file => $seqfile, 
								-format => $format);
my $seq_obj = $seqio_obj->next_seq;
my @pro_seqs = ();

push @pro_seqs, $seq_obj->translate(-frame => 0)->seq;
push @pro_seqs, $seq_obj->translate(-frame => 1)->seq;
push @pro_seqs, $seq_obj->translate(-frame => 2)->seq;
push @pro_seqs, $seq_obj->revcom->translate(-frame => 0)->seq;
push @pro_seqs, $seq_obj->revcom->translate(-frame => 1)->seq;
push @pro_seqs, $seq_obj->revcom->translate(-frame => 2)->seq;

my $protein = 0;
my $prolen = 0;
foreach my $pro_seq (@pro_seqs){
	my @orfs = split(/\*/, $pro_seq);
	foreach my $orf (@orfs){
		if(length($orf) > $prolen){
			$protein = $orf;
			$prolen = length($orf);
		}
	}
}

$protein =~ s/^[A-KN-Z]+M/M/;
$protein = $protein."*";
$prolen = length($protein);

print "$prolen\n$protein\n";