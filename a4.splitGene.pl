###########################################################
###	@Author		Xi Gu									###
###	@Time		Oct 28, 2014							###
###	@Function	Split Query into Multiple Files			###
###########################################################

#!/usr/bin/perl
# run the script: time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/01.Protein/run.0 1000

use strict;

my $seqFile = shift @ARGV;			# proteome seq file
my $refFile = shift @ARGV;			# gene id file
my $outFolder = shift @ARGV;			# output file
my $splitSize = shift @ARGV;


open(SEQ, $seqFile) or die "Cannot open $seqFile: $!";

my %hash = ();
my $gene = 0;
while (my $line = <SEQ>){
	chomp $line;
	if($line =~ /^>/){
		$gene = $line;
		$gene =~ s/>//;
		$hash{$gene} = ();
	}else{
		push @{$hash{$gene}},$line;
	}
}

open(REF, $refFile) or die "Cannot open $refFile: $!";
my $i = 1;
my $group = $i * $splitSize;
system("mkdir -p $outFolder/$group");
open(TGT, ">$outFolder/$group/$group.fasta")or die "Cannot write $group.fasta: $!";

foreach my $line (<REF>){
	chomp $line;
	my @lines = split(/\t/, $line);
	$line = $lines[0];
	if($line eq ""){
		$i = $i + 1;
		close TGT;
		$group = $i * $splitSize;
		system("mkdir -p $outFolder/$group");
		open(TGT, ">$outFolder/$group/$group.fasta")or die "Cannot write $group.fasta: $!";
		next;
	}
	print TGT ">$line\n";
	print TGT join("\n", @{$hash{$line}}), "\n";
}

if(-z "$outFolder/$group/$group.fasta"){
	system("rm -r $outFolder/$group");
}

close REF;
close SEQ;
close TGT;