#!/usr/bin/perl -w
# run the script: time perl 00.script/d17.transform.ID.pl 01.data/06.TargetTranscriptome/StHe.transcriptome.v0.fa 01.data/06.TargetTranscriptome/StHe.transcriptome.v1.fa

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
		$line =~ s/^>//;
		my @lines = split(/\s+/, $line);
		my $ID = shift @lines;
		my $last = pop @lines;
		if($last =~ /^[0-9]+$/){
			$ID = $ID.".$count";
		}
		$line = ">".$ID." ".join(" ", @lines)." ".$last;
	}
	print TGT "$line\n";
}

close SRC;
close TGT;