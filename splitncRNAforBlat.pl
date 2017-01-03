#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $outfolder = shift @ARGV;
my $out_prefix = shift @ARGV || $srcfile;
my $size = shift @ARGV;

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => 'fasta');

my $i = 1;
my $count = 1;
open(my $io, ">$outfolder/$out_prefix.$i.fasta");
print "Writing part $i...\n";

while(my $seqio = $seqio_obj->next_seq){
	$count ++;
	if($count % $size == 0){
		close($io);
		$i ++;
		print "Writing part $i...\n";
		open($io, ">$outfolder/$out_prefix.$i.fasta");
	}
	print $io ">", $seqio->id, " ", $seqio->description, "\n", $seqio->seq, "\n";
}

close($io);