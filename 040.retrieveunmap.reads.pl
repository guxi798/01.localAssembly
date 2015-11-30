#!/usr/bin/perl -w
# run the script: time perl 00.script/040.retrieveunmap.reads.pl 03.blast/03.bowtie.nucl/run.0 B1 R1 01.data/02.Fasta/B1/B1.R1.fasta_simple.fasta >> 04.retrieve.reads/03.bowtie.nucl/run.0/99999/retrieved.99999.R1.fasta
# under master script 04.folder.retrievebowtie.reads.pl

use strict;

my $srcfolder = shift @ARGV;
my $sam = shift @ARGV;
my $mode = shift @ARGV;
my $reffile = shift @ARGV;

opendir(SRC, "$srcfolder") or die "ERROR: Cannot open $srcfolder: $!";
open(REF, "$reffile") or die "ERROR: Cannot open $reffile: $!";

my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));

my %hash = ();
foreach my $sub (@subs){
	my $file = 0;
	if($mode eq "R1"){
		$file = "$srcfolder/$sub/bowtie.out.$sub.$sam.R1.tab";
	}elsif($mode eq "R2"){
		$file = "$srcfolder/$sub/bowtie.out.$sub.$sam.R2.tab";
	}elsif($mode eq "S"){
		$file = "$srcfolder/$sub/bowtie.out.$sub.$sam.single.tab";
	}elsif($mode eq "F"){
		$file = "$srcfolder/$sub/bowtie.out.$sub.$sam.tab";
	}
	open(SUB, $file) or die "ERROR: Cannot open $file: $!";
	
	foreach my $line (<SUB>){
		chomp $line;
		my @lines = split(/\t/, $line);
		if(not exists $hash{$lines[0]}){
			$hash{$lines[0]} = 0;
		}
	}
	
	close SUB;
}

my $prev = 0;
my $mark = 0;
my @keys = sort {$a <=> $b} keys %hash;
my $key = shift @keys;

while(my $line = <REF>){
	chomp $line;
	if($line =~ />/){
		$line =~ s/>//;
		if(! defined $key){
			$mark = 1;
			print ">$line\n";
			next;
		}
		if($line > $prev and $line < $key){
			$mark = 1;
			print ">$sam","_","$line\n";
		}elsif($line > $key){
			for($prev > $key){
				$prev = $key;
				$key = shift @keys;
			}
			if(defined $key and $line > $prev and $line < $key or ! defined $key){
				$mark = 1;
				print ">$sam","_","$line\n";
			}
		}else{
			$mark = 0;
		}
	}else{
		if($mark){
			print "$line\n";
		}
	}
}

close REF;