#!/usr/bin/perl -w
# run the script: time perl 00.script/b2.extract.gene.of.interest.pl 
#					01.data/00.PriorData/SUT.txt 
#					07.map.back/03.bowtie.nucl 
#					06.assembly/03.bowtie.nucl 
#					01.data/04.GeneOfInterest/GeneID.v1.txt
#					11.stable.run/SUT.blast.summary.txt
#					11.stable.run/SUT.contig.seq.fasta

use strict;

my $reffile = shift @ARGV;			# file contains genes of interests. e.g. SUT
my $blastfolder = shift @ARGV;		# blast result folder
my $trinityfolder = shift @ARGV;	# trinity result folder
my $infofile = shift @ARGV;			# gene info file contains protein length and group ID
my $tgtfile1 = shift @ARGV; 		# overall summary output
my $tgtfile2 = shift @ARGV; 		# include detailed sequences

my @path = split(/\//, $tgtfile1);
system("mkdir -p $path[0]");

open(REF, $reffile) or die "ERROR: Cannot open $reffile: $!";
opendir(SRC, $blastfolder) or die "ERROR: Cannot open $blastfolder: $!";
open(INF, $infofile) or die "ERROR: Cannot open $infofile: $!";

#### based on info file, construct a hash reference
my %hash = ();
foreach my $line (<INF>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}){
		$hash{$lines[0]} = [$lines[1], $lines[5], $lines[6]]; # protein length, group ID, expression
	}
}

#### get how many runs
my @runs = grep(/run\.[0-9]+$/, readdir(SRC));
@runs = map {$_ =~ s/run\.//g; $_} @runs;
@runs = sort {$a <=> $b} @runs;
#my @runs = ("run.9");
#### get blast result records for each gene in each run
system("> $tgtfile1");
foreach my $line (<REF>){	# for each GOI
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}){
		next;
	}
	my $prolen = ${$hash{$lines[0]}}[0];
	my $group = ${$hash{$lines[0]}}[1];
	#my $expr = ${$hash{$lines[0]}}[2];
	#system("echo \">$lines[0]	$lines[1]	ProLen=$prolen	$group $expr\" >> $tgtfile1");
	system("echo \">$lines[0]	$lines[1]	ProLen=$prolen	$group\" >> $tgtfile1");
	foreach my $run (@runs){	# for each run
		system("echo \"\@run.$run\" >> $tgtfile1");
		system("grep \"$lines[0]\" $blastfolder/run.$run/$group/$group.contigs.blast.out >> $tgtfile1");
	}
}

#### get contig name from the output file from previous step
open(SUM, $tgtfile1) or die "ERROR: Cannot open $tgtfile1: $!";
#### open sequence output file
open(TGT, ">$tgtfile2") or die "ERROR: Cannot write $tgtfile2: $!";

my $group = 0;		# group ID
my $run = 0;		# run ID
my $contig = 0;		# contig ID
my %contig = ();	# hash to store contig seqs
foreach my $line (<SUM>){	# for each blast record, corresponding to each contig
	chomp $line;
	my @lines = split(/\s+/, $line);
	if($line =~ /^>/){		# if it is the info line about GOI
		$group = $lines[3];	# get group ID to open corresponding file
		print TGT "$line\n";
	}elsif($line =~ /^@/){	# if it is the info line about run ID
		$run = $line;
		$run =~ s/@//;
		print TGT "$run\n";
		$contig = 0;		# initialize contig structure to store contig seqs
		%contig = ();
		# open contig file, that is the trinity output file, group ID and run ID are needed.
		open(SEQ, "$trinityfolder/$run/$group/Trinity.new.fasta") or die "ERROR: Cannot open $trinityfolder/$run/$group/Trinity.new.fasta: $!";
		foreach my $row (<SEQ>){
			chomp $row;
			if($row =~ /^>/){
				my @rows = split(/\s+/, $row);
				$contig = $rows[0];
				$contig =~ s/^>//;
			}else{
				if(not exists $contig{$contig}){
					$contig{$contig} = [$row];
				}else{
					push @{$contig{$contig}}, $row;
				}
			}
		}
		close SEQ;
	}else{					# if it is the blast record line
		my @ids = split(/\|/, $lines[0]);	# get contig ID
		my $id = $ids[0];
		if(exists $contig{$id}){
			print TGT ">$id\n";
			print TGT join("\n", @{$contig{$id}}), "\n";
		}		
	}
}

close REF;
close SRC;
close INF;
close TGT;
