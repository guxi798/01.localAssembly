#!/usr/bin/perl
# run the script: time perl 00.script/06.truncate.header.pl 06.assembly/03.bowtie.nucl/run.0/33000/Trinity.old.fasta 06.assembly/03.bowtie.nucl/run.0/33000/Trinity.fasta

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile) or die "Error: Cannot open $srcfile: $!";
open(TGT, ">$tgtfile") or die "Error: Cannot open $tgtfile: $!";

my @seq = ();
my $contig = 0;
my $pre = 0;
my $gc = 0;
my $len = 0;
foreach my $line (<SRC>){
	chomp $line;
	if($line =~ /^>/){
		my @lines = split(/\s+/, $line);
		if($contig ne 0){
			my $pergc = sprintf "%.2f", $gc/$len*100;
			$pre = $contig." GC=".$pergc;
			print TGT "$pre\n";
			print TGT join("\n", @seq), "\n";
		}
		$contig = $lines[0]." ".$lines[1];
		@seq = ();
		$gc = 0;
		$len = 0;
	}else{
		push @seq, $line;
		$gc = $gc + ($line =~ tr/gcGC//);
		$len = $len + length($line);
	}
}

my $pergc = sprintf "%.2f", $gc/$len*100;
$pre = $contig." GC=".$pergc;
print TGT "$pre\n";
print TGT join("\n", @seq), "\n";

close SRC;
close TGT;