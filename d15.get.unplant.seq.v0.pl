#!/usr/bin/perl -w
# run the script: time perl ../../00.script/d15.get.unplant.seq.pl tblastn.batch transcriptome.part3.v3.fa ../04.unplant/transcriptome.part3.vNA.fa

use strict;
use Bio::SeqIO;

my $reffolder = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

opendir(REF, $reffolder);
my @files = sort(grep(/^rh/, readdir(REF)));
closedir REF;

my %hash = ();
foreach my $file (@files){
	open(FHL, "$reffolder/$file/$file.out");
	while(my $line = <FHL>){
		chomp $line;
		my @lines = split(/\s+/, $line);
		if(not exists $hash{$lines[1]}){
			$hash{$lines[1]} = 0;
		}
	}
	close FHL;
}

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");
open(TGT, ">$tgtfile");
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
	}
}
close TGT;

print "There are $count1 contigs in total\n";
print "$count2 contigs are not mapped to plants, and will be blast against fungi and bacteria RNA\n";