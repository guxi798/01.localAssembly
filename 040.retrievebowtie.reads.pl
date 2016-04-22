#!/usr/bin/perl -w
# run the script: time perl 00.script/040.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.0/1000/bowtie.out.1000.X88.R1.tab 01.data/02.Fasta/X88/X88.R1.fasta_simple.fasta 1000 >> 04.retrieve.reads/03.bowtie.nucl/run.0/1000/retrieved.1000.R1.fasta
# under master script 04.folder.retrievebowtie.reads.pl

use strict;

my $srcfile = shift @ARGV;
my $reffile = shift @ARGV;
my $blocksize = shift @ARGV;

my @path1 = split(/\//, $reffile);
my @path2 = split(/\./, pop @path1);
my $sam = shift @path2;

open(SRC, "$srcfile") or die "ERROR: Cannot open $srcfile: $!";
open(REF, "$reffile") or die "ERROR: Cannot open $reffile: $!";

my $run = 0;
my $previous = 0;
my $seq = 0;
my $count = 0;
my $refline = 0;
my $line = 0;

while (1){
	my @records = ();
	my @ids = ();
	$run ++;
	#print "$run\n";
	
	while($refline = <REF>){
		if($refline =~ /^>/){
			$count ++;
			chomp $refline;
			$refline =~ s/^>//;
			$refline =~ s/\/1$//;
			$refline =~ s/\/2$//;
			my @reflines = split(/\s+/, $refline);
			push @ids, $reflines[0];
			if(!$seq){next;}
			push @records, $seq;
			$seq = 0;
			next;
		}
		chomp $refline;
		if($seq){
			$seq = $seq.$refline;
		}else{
			$seq = $refline;
		}
		
		if($count >= $blocksize){
			$count = 0;
			last;
		};
	}
	
	if($refline){
		push @records, $seq;
		$seq = 0;
	}

	while(1){
		if($previous != 0){
			if($previous > ($run - 1) * $blocksize and $previous <= $run * $blocksize){
				#print ">$sam","_","$previous\n";
				print ">$ids[$previous % $blocksize - 1]\n";
				print "$records[$previous % $blocksize - 1]\n";
				$previous = 0;
			}else{
				last;
			}
		}
		
		$line = <SRC>;
		if(!$line){last;}
		
		chomp $line;
		my @lines = split(/\t/, $line);
		my $pos = $lines[0] % $blocksize;
		if($pos == 0){$pos = $blocksize;}
		if($lines[0] > $run * $blocksize){
			$previous = $lines[0];
			last;
		}else{
			#print ">$sam","_","$lines[0]\n";
			print ">$ids[$pos - 1]\n";
			print "$records[$pos - 1]\n";
		}
	}
	
	if(!$line){
		last;
	}

}

close SRC;
close REF;