#!/usr/bin/perl -w
# run the script: time perl 00.script/a7.retrieve.dustMasker.fasta.pl 01.data/01.Fastq/OrAe1G/OrAe1G.R1.fastq 01.data/02.Fasta/OrAe1G/OrAe1G.R1.masked.fasta 01.data/01.Fastq/OrAe1G/OrAe1G.R1.clean.fastq 0.5
use strict;

my $reffile = shift @ARGV;
my $mskfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $cut = shift @ARGV;

open(REF, $reffile);
open(MSK, $mskfile);
open(TGT, ">$tgtfile");
open(TMP, ">$tgtfile.tmp");

while(my $refline = <REF>){
	chomp $refline;
	if($refline =~ /^@/){
		my $id = <MSK>;
		chomp $id;
		my $seq = <MSK>;
		chomp $seq;
		my $test = 0;
		while($test = <MSK>){
			if($test =~ /^>/){last;}
			chomp $test;
			$seq = $seq.$test;
		}
		if($test){
			seek(MSK, -length($test), 1);
		}
		
		my $comp1 = $refline;
		my $comp2 = $id;
		$comp1 =~ s/^@//;
		$comp2 =~ s/^>//;
		if($comp1 ne $comp2){
			print "something went wrong: $comp1\t$comp2\n";
			die;
		}
		
		my $len = length($seq);
		my $mask = $seq =~ tr/a-z//;
		if($mask/$len < $cut){
			print TGT "$refline\n";
			$refline = <REF>;
			print TGT "$refline";
			$refline = <REF>;
			print TGT "$refline";
			$refline = <REF>;
			print TGT "$refline";
		}else{
			print TMP "$id\n";
			print TMP "$seq\n";
			<REF>;
			<REF>;
			<REF>;
		}
	}
}

close REF;
close MSK;
close TGT;
close TMP;
