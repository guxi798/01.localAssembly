#!/usr/bin/perl -w
# run the script: time perl ../../00.script/d16.remove.nonplant.seq.pl blastn.batch transcriptome.part3.v3.fa ../../01.data/06.TargetTranscriptome/transcriptome.part4.v4.fa transcriptome.part4.vNA.fa

use strict;
use Bio::SeqIO;

my $reffolder = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $tgtfile2 = shift @ARGV;

opendir(REF, $reffolder);
my @files = sort(grep(/^rh/, readdir(REF)));
closedir REF;

my %hash = ();
foreach my $file (@files){
	print "Processing $file...\n";
	open(SRC, "$reffolder/$file/$file.hash.out");
	foreach my $line(<SRC>){
		chomp $line;
		print "$line\n";
		my @lines = split(/\s+/, $line);
		if(not exists $hash{$lines[0]}){
			$hash{$lines[0]} = 0;
		}
	}
}

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");
open(TGT, ">$tgtfile");
open(TGT2, ">$tgtfile2");
my $count1 = 0;
my $count2 = 0;
while(my $seq_obj = $seqio_obj->next_seq){
	$count1 ++;
	my $line = $seq_obj->id;
	my $desc = $seq_obj->desc;
	my $seq = $seq_obj->seq;
	
	if(not exists $hash{$line}){
		$count2 ++;
		print TGT ">$line $desc\n";
		print TGT "$seq\n";
	}else{
		print TGT2 ">$line $desc\n";
		print TGT2 "$seq\n";
	}
}
close TGT;
close TGT2;

print "There are $count1 contigs in total\n";
print $count1-$count2, " contigs are mapped to bacteria or fungi\n";
print "$count2 contigs are retained as plant sequences\n";

## function
sub min{
	my $l = shift @_;
	my $r = shift @_;
	if($l <= $r*0.7){return $r;}
	if($l <= $r){return $l}
	else{return $r;}
}