#!/usr/bin/perl -w
# run the script: time perl 00.script/13.prepare.file.pl 01.data/05.splitGenes/02.Transcript/run.3/contigs.good.fasta 01.data/05.splitGenes/02.Transcript/run.3/contig2gene.txt

use strict;

my $reffile = shift @ARGV;	
my $outfile = shift @ARGV;
my $genefile = shift @ARGV;		

open(SRC, $reffile) or die "ERROR: Cannot open $reffile: $!";
open(TGT1, ">$outfile");
open(TGT2, ">$genefile");

my $count = 0;
my $pregene = '';

foreach my $line (<SRC>){	
	chomp $line;
	if($line =~ /^>/){		
		$line =~ s/>//;
		my @lines = split(/ /, $line);
		my @part = split(/_/, $lines[0]);
		my $gene = join("_", @part[0..1]);
		my $contig = 0;
		
		if($gene eq $pregene){
			$contig = "c".$count."_$part[1]_$part[2]";
		}else{
			$count ++;
			$contig = "c".$count."_$part[1]_$part[2]";
		}
		print TGT1 ">$contig $lines[1]\n";
		print TGT2 "$contig\t$lines[1]\n";
		
		$pregene = $gene;
	}else{
		print TGT1 "$line\n";
	}

}

close SRC;
close TGT1;
close TGT2;
