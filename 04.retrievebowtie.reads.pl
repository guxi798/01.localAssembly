#!/usr/bin/perl -w
# run the script: time perl 00.script/04.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl 1000 tgtfolder
# under master script 04.folder.retrievebowtie.reads.pl

use strict;

my $srcfile = shift @ARGV;
my $mode = shift @ARGV;
my $tgtfile1 = shift @ARGV;

open(SRC, "$srcfile");
open(TGT1, ">$tgtfile1");
if(!($srcfile =~ /long\.sam$/)){
	my $tgtfile2 = shift @ARGV;
	open(TGT2, ">$tgtfile2");
}

if($mode eq "both"){
	my $unmapfile1 = shift @ARGV; 
	open(UNM1, ">$unmapfile1");
	
	if(!($srcfile =~ /long\.sam$/)){
		my $unmapfile2 = shift @ARGV;
		open(UNM2, ">$unmapfile2");
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
	
	if($srcfile =~ /long\.sam$/){
		if($decnum == 0){
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
		if(($decnum == 0) or ($hexflag =~ /[4|5|6|7][1|3|9|b]$/)){
			print TGT1 ">$readid/1\n";
			print TGT1 "$seq\n";
		}elsif($hexflag =~ /[8|9|a|b][1|3|9|b]$/){
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
if(!($srcfile =~ /long\.sam$/)){
	close TGT2;
}
close UNM1;
close UNM2;
#close REC1;
#close REC2;


	
