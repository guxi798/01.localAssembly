#!/usr/bin/perl -w
# run the script: time perl 00.script/z2.getseq.with.geneID.pl 01.data/00.PriorData/tubulin.txt 01.data/00.PriorData/transcriptome.fa 01.data/05.splitGenes/02.Transcript/run.0/1000/1000.fasta

use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $seqfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile);
my %hash = ();
while(my $line = <SRC>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[1]}){
		$hash{$lines[1]} = 0;
	}
}
close SRC;

open(TGT, ">$tgtfile");
my $seqio_obj = Bio::SeqIO->new(-file => $seqfile, -format => "fasta");
while(my $seqobj = $seqio_obj->next_seq){
	my $id = $seqobj->id;
	if(exists $hash{$id}){
		my $seq = $seqobj->seq;
		print TGT ">$id\n";
		print TGT "$seq\n";
	}
}

close TGT;