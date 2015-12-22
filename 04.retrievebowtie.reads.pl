#!/usr/bin/perl -w
# run the script: time perl 00.script/04.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl 1000 tgtfolder
# under master script 04.folder.retrievebowtie.reads.pl

use strict;

my $srcfile = shift @ARGV;
my $mode = shift @ARGV;
my $tgtfile1 = shift @ARGV;
my $unmapfile1 = shift @ARGV; 
my $tgtfile2 = shift @ARGV;
my $unmapfile2 = shift @ARGV;
my ($run) = $srcfile =~ /run\.([0-9]+)/;
my $oldrun = $run - 1;

if($run != 0){
	my $oldfile1 = $tgtfile1;
	$oldfile1 =~ s/run\.$run/run\.$oldrun/;
	open(OLD1, $oldfile1);
}

open(SRC, "$srcfile");
open(TGT1, ">$tgtfile1");
if($tgtfile2){
	open(TGT2, ">$tgtfile2");
	if($run != 0){
		my $oldfile2 = $tgtfile2;
		$oldfile2 =~ s/run\.$run/run\.$oldrun/;
		open(OLD2, $oldfile2);
	}
}

if($mode eq "both"){
	open(UNM1, ">$unmapfile1");
	
	if($unmapfile2){
		open(UNM2, ">$unmapfile2");
	}
}

## if not the first run, then store previous reads into hash, and print to result file
## if current run recover old reads, ignore, if find new reads, add to current result file
## this way will ensure old read information would not be lost
my %hash1 = ();
my %hash2 = ();
if($run != 0){
	while(my $line = <OLD1>){
		chomp $line;
		$line =~ s/^>//;
		$line =~ s/\/1$//;
		my $seq = <OLD1>;
		if(not exists $hash1{$line}){
			chomp $seq;
			$hash1{$line} = $seq;
			print TGT1 ">$line/1\n$seq\n";
		}
	}
	
	if($tgtfile2){
		while(my $line = <OLD2>){
			chomp $line;
			$line =~ s/^>//;
			$line =~ s/\/2$//;
			my $seq = <OLD2>;
			if(not exists $hash2{$line}){
				chomp $seq;
				$hash2{$line} = $seq;
				print TGT2 ">$line/2\n$seq\n";
			}
		}
	}
}

while (my $line = <SRC>){
	if($line =~ /^@/) {next;}
	chomp $line;
	
	my @lines = split(/\t/, $line);
	my $readid = $lines[0];
	my $flag = $lines[1];
	my $seq = $lines[9];
	my $hexnum = sprintf("0x%x", $flag);
	my $decnum = sprintf("%d", $flag);
	my $hexflag = "$hexnum";
	
	if(!$tgtfile2){
		if($decnum == 0 and not exists $hash1{$readid}){
			print TGT1 ">$readid/1\n";
			print TGT1 "$seq\n";
		}else{
			if($mode ne "both"){next;}
			if($hexflag =~ /4$/){
				print UNM1 ">$readid/1\n";
				print UNM1 "$seq\n";
			}else{
				print "ERROR: the format is supposed to be single!\n";
			}
		}
	}else{
		if(($decnum == 0 and not exists $hash1{$readid}) or ($hexflag =~ /[4|5|6|7][1|3|9|b]$/ and not exists $hash1{$readid})){
			print TGT1 ">$readid/1\n";
			print TGT1 "$seq\n";
		}elsif($hexflag =~ /[8|9|a|b][1|3|9|b]$/ and not exists $hash2{$readid}){
			print TGT2 ">$readid/2\n";
			print TGT2 "$seq\n";
		}else{
			if($mode ne "both"){next;}
			if(($hexflag =~ /[4|5|6|7][5|d]$/) or ($hexflag =~ /4$/)){
				print UNM1 ">$readid/1\n";
				print UNM1 "$seq\n";
			}elsif($hexflag =~ /[8|9|a|b][5|d]$/){
				print UNM2 ">$readid/2\n";
				print UNM2 "$seq\n";
			}
		}
	}
}

close SRC;
close TGT1;
if($tgtfile2){close TGT2;}
if($unmapfile1){close UNM1;}
if($unmapfile2){close UNM2;}


	
