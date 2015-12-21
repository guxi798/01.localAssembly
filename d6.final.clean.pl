#!/usr/bin/perl -w
# run the script: perl 00.script/d6.final.clean.pl ../../12.Parasite/01.OrAe/01.data/01.Fastq/OrAe2G/OrAe2G.R1.fastq 01.data/01.Fastq/OrAe2G/OrAe2G.R1.fastq 01.data/01.Fastq/OrAe2G/OrAe2G.R1.clean

use strict;

my $reffile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(REF, $reffile);
open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my %hash = ();
foreach my $line(<REF>){
	if(!($line =~ /^\@HW|^\@SO|^\@IL/)){next;}
	chomp $line;
	if(not exists $hash{$line}){
		$hash{$line} = 0;
	}
}

while (my $line = <SRC>){
	if($line =~ /^\@HW|^\@SO|^\@IL/){
		chomp $line;
		if(exists $hash{$line}){
			print TGT "$line\n";
			$line = <SRC>;
			print TGT "$line";
			$line = <SRC>;
			print TGT "$line";
			$line = <SRC>;
			print TGT "$line";
		}
	}
}

close REF;
close SRC;
close TGT;