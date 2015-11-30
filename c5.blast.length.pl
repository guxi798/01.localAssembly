#!/usr/bin/perl -w
# run the script: time perl 00.script/c3.combine.fasta.pl 06.assembly/03.bowtie.nucl/run.0

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;

## start running the script
opendir(SRC, $srcfolder) or die "Cannot open $srcfolder: $!";
my @subs = sort(grep(/[0-9]+/, readdir(SRC)));
my $out = "$srcfolder/length.summary.txt";
open(TGT, ">$out");
print TGT "blast_contig_len\tblast_ref_len\n";

my $contig_len = 0;
my $ref_len = 0;
my $count = 0;
my $evalue = 0;
my $mismatch = 0;
foreach my $sub (@subs){
	open(SUB, "$srcfolder/$sub/$sub.contigs.blast.out");

	foreach my $line (<SUB>){
		$count ++;
		chomp $line;
		my @lines = split(/\t/, $line);
		print TGT abs($lines[7] - $lines[6]), "\t", abs($lines[9] - $lines[8]), "\n";
		$contig_len += abs($lines[7] - $lines[6]);
		$ref_len += abs($lines[9] - $lines[8]);
		$evalue += $lines[10];
		$mismatch += $lines[4];
	}
	close SUB;
}
$contig_len = $contig_len/$count;
$ref_len = $ref_len/$count;
$evalue = $evalue/$count;
$mismatch = $mismatch/$count;
print "$contig_len\t$ref_len\t$evalue\t$mismatch\n";

close SRC;
close TGT;
