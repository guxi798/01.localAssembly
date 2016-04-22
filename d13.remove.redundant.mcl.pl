#!/usr/bin/perl -w
# run the script: time perl 00.script/d13.remove.redundant.mcl.pl 12.MCL.remove.redundancy/mcl.out.txt 12.MCL.remove.redundancy/transcriptome.part3.v2.fa 01.data/06.TargetTranscriptome/transcriptome.part3.v3.fa

use strict;
use Bio::SeqIO;

my $reffile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");

my %hash = ();
while(my $seq_obj = $seqio_obj->next_seq){
	my $line = $seq_obj->id;
	my @temp = split(/\s+/, $seq_obj->desc);
	my $len = shift @temp;
	$len =~ s/len=//;
	my $seq = $seq_obj->seq;
	
	if(not exists $hash{$line}){
		$hash{$line} = [$len, $seq, \@temp];
	}else{
		die "Error: $line occurs more than once\n";
	}
}

open(REF, $reffile);
open(TGT, ">$tgtfile");

my $count = 0;
foreach my $line (<REF>){
	my $longestlen = 0;
	my $longestID = 0;
	my $longestseq = 0;
	my $longesttemp = 0;

	$count++;
	chomp $line;
	my @lines = split(/\s+/, $line);
	foreach my $ID (@lines){
		if(exists $hash{$ID}){
			my $len = ${$hash{$ID}}[0];
			if($len > $longestlen){
				$longestlen = $len;
				$longestID = $ID;
				$longestseq = ${$hash{$ID}}[1];
				$longesttemp = ${$hash{$ID}}[2];
			}
		}else{
			die "Error: $ID doesn't exist in the sequence file.\n";
		}
	}
	 
	print "Group_$count\t$longestID\tlen=$longestlen\n";
	print TGT ">$longestID len=$longestlen ", join(" ", @{$longesttemp}), "\n";
	print TGT "$longestseq\n";
}

close REF;
close TGT;

