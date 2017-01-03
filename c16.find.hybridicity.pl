#!/usr/bin/perl

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");
my $i = 0;
my $j = 0;
my $mismatch_l = 0;
my $mismatch_t = 0;
my $bits_l = 0;
my $bits_t = 0;
my $gene = 0;

while(my $line = <SRC>){
	chomp $line;
	if($line =~ /^>/){
		if($i == 0 and $j == 0){next;}
		
		#print TGT "$gene\t", $mismatch_l/$i, "\t", $mismatch_t/$j, "\t", $bits_l/$i, "\t", $bits_t/$j, "\n";
		#print TGT "$gene\t$i\t$mismatch_l\t$bits_l\t$j\t$mismatch_t\t$bits_t\n";
		print TGT "$gene\t$mismatch_l\t$bits_l\t$mismatch_t\t$bits_t\n";
		$mismatch_l = $mismatch_t = $bits_l = $bits_t = $i = $j = 0;
	}else{
		my @lines = split(/\s+/, $line);
		$gene = $lines[2];
		if($lines[0] eq "Local" and $lines[12] > $bits_l){
			$i ++;
			#$mismatch_l += $lines[5];
			#$bits_l += $lines[12];
			$mismatch_l = $lines[5];
			$bits_l = $lines[12];
		}elsif($lines[0] eq "Trinity" and $lines[12] > $bits_t){
			$j ++;
			#$mismatch_t += $lines[5];
			#$bits_t += $lines[12];
			$mismatch_t = $lines[5];
			$bits_t = $lines[12];
		}
	}
}

close SRC;
close TGT;