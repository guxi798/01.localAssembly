#!/usr/bin/perl -w
# run the script: time perl 00.script/z1.transform.ID.pl 01.data/00.PriorData/total.RNA.part5.fasta 01.data/00.PriorData/total.RNA.part51.fasta

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

my $count = 0;
open(SRC, $srcfile);
open(TGT, ">$tgtfile");

while(my $line = <SRC>){
	chomp $line;
	if($line =~ /^>/){
		$count ++;
		$line = '>'.$count;
	}
	print TGT "$line\n";
}

close SRC;
close TGT;