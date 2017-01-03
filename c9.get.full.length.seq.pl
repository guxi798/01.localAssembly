#!/usr/bin/perl -w
## run the script: time perl 00.script/c9.get.full.length.seq.pl 08.full.length/count3 08.full.length/full.length.contigs.fasta 08.full.length/Final.fasta

use strict;

my $reffile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(REF, $reffile);
my %hash = ();
foreach my $line (<REF>){
	chomp $line;
	my @lines = split(/\t/, $line);
	my $gene = shift @lines;
	my $contig = shift @lines;
	my $run = shift @lines;
	if(not exists $hash{$run}{$gene}{$contig}){
		$hash{$run}{$gene}{$contig} = 0;
	}else{
		die "Something unexpected occurs.\n";
	}
}
close REF;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my $count = 0;
while(my $line = <SRC>){
	chomp $line;
	if($line =~ /^>/){
		$line =~ s/^>//;
		my @lines = split(/\s+/, $line);
		my $contig = $lines[0];
		my $seqlen = $lines[1];
		my $gene = $lines[2];
		my $run = pop @lines;
		$run = $run - 1;
		#print "$run\t$group\t$contig\n";
		if(exists $hash{$run}{$gene}{$contig}){
			$count ++;
			#print "I'm here\n";
			print TGT ">$contig.$count $seqlen $gene $run\n";
			$line = <SRC>;
			print TGT "$line";
		}
	}
}

close SRC;
close TGT;