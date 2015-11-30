#!/usr/bin/perl
# run the script: time perl 00.script/b1.change.geneID.pl 01.data/00.PriorData/Ptrichocarpa_210_v3.0.protein_primaryTranscriptOnly.fa 01.data/00.PriorData/proteome.fa

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

foreach my $line (<SRC>){
	chomp $line;
	if($line =~ /^>/){
		$line =~ s/>//;
		my @lines = split(/\s+/, $line);
		my @ids = split(/\./, $lines[0]);
		pop @ids;			# remove transcript index
		my $id = join(".", @ids);
		print TGT ">$id\n";
	}else{
		print TGT "$line\n";
	}
}

close SRC;
close TGT;