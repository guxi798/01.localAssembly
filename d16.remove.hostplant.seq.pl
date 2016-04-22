#!/usr/bin/perl -w
# run the script: time perl ../../00.script/d16.remove.hostplant.seq.pl transcriptome.part1.v2.blastn.xml.out transcriptome.part1.v2.fa ../../01.data/06.TargetTranscriptome/transcriptome.part1.v4.fa transcriptome.part1.vNA.fa 0.8

use strict;
use Bio::SeqIO;
use Bio::SearchIO;

my $reffile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $tgtfile2 = shift @ARGV;
my $cut = shift @ARGV;

my %hash = ();
my $blast = new Bio::SearchIO(-format => 'blastxml', -file => "$reffile");
while(my $result = $blast->next_result){
	my @names = split(/\s+/, $result->query_description);
	my $qname = $names[0];
	my $qlen = $result->query_length;
	while(my $hit = $result->next_hit){
		#my $frac_iden = $hit->frac_identical;
		#if($frac_iden < 0.8){next;}
		my $hitname = $hit->name;
		my $hitlen = $hit->length;
		my $len = &min($qlen, $hitlen);
		#my $len = $hitlen;
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		my $frac = $identity/$len;
		#print "$hitname\t$qlen\t$hitlen\t$identity\t$frac\n";
		if($frac > $cut){
			if(not exists $hash{$qname}){
				$hash{$qname} = 0;
				print "$hitname\t$qlen\t$hitlen\t$identity\t$frac\n";
			}
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