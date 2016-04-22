#!/usr/bin/perl -w
# after running 10.transfer.sature.seq.pl
# run the script: time perl 00.script/d11.remove.redundant.contigs.pl  01.data/06.TargetTranscriptome/transcriptome.part1.v1.fa  01.data/06.TargetTranscriptome/transcriptome.part1.v2.fa

use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my $prevlen = 0;
my $prevID = "";
my $prevline = 0;
my $prevseq = "";
my $count1 = 0;
my $count2 = 0;

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");

while(my $seq_obj = $seqio_obj->next_seq){
	$count1 ++;
	my $line = $seq_obj->id;
	my @temp = split(/\s+/, $seq_obj->desc);
	my @ID = split(/_/, $line);
	pop @ID;
	my $ID = join("_", @ID);
	my $len = $temp[0];
	my $run = pop @temp;
	$len =~ s/len=//;
	my $seq = $seq_obj->seq;
	
	if($ID ne $prevID and $prevID ne ""){
		$count2 ++;
		print TGT ">$prevline\n";
		print TGT "$prevseq\n";
	}else{
		if($len <= $prevlen){
			next;
		}
	}
	
	$prevlen = $len;
	$prevID = $ID;
	$prevline = $line." ".$seq_obj->desc;
	$prevseq = $seq;
}

$count2 ++;
print TGT ">$prevline\n";
print TGT "$prevseq\n";

print "There are $count1 contigs before removing redundancy\n";
print "There are $count2 contigs after removing redundancy\n";

close SRC;
close TGT;