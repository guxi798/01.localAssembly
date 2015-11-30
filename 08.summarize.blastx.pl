#!/usr/bin/perl -w
# run the script: time perl 00.script/08.summarize.blastx.pl 
# 07.map.back/03.bowtie.nucl/run.0/2000/2000.contigs.blast.out 
# 07.map.back/03.bowtie.nucl/run.1/2000/2000.contigs.blast.out 
# 08.evaluate/03.bowtie.nucl/run.1/2000/2000.compare.tab

use strict;

my $oldfile = shift @ARGV;
my $newfile = shift @ARGV;
my $outfile = shift @ARGV;

open(OLD, $oldfile);
open(NEW, $newfile);
open(TGT, ">$outfile");

my %hash = ();

foreach my $line (<OLD>){
	chomp $line;
	my @lines = split(/\t/, $line);
	my $gene = $lines[1];
	my $identity = $lines[2];
	my $len = $lines[3];
	my $evalue = $lines[10];
	if(not exists $hash{$gene}){
		$hash{$gene} = [$identity, $len, $evalue, "NA", "NA", "NA"];
	}else{
		if((${$hash{$gene}}[1] < $len) or (${$hash{$gene}}[1]==$len and ${$hash{$gene}}[0]<$identity)){
			$hash{$gene} = [$identity, $len, $evalue, "NA", "NA", "NA"];
		}
	}
}

foreach my $line (<NEW>){
	chomp $line;
	my @lines = split(/\t/, $line);
	my $gene = $lines[1];
	my $identity = $lines[2];
	my $len = $lines[3];
	my $evalue = $lines[10];
	if(exists $hash{$gene}){
		if(${$hash{$gene}}[3] eq "NA"){
			${$hash{$gene}}[3] = $identity;
			${$hash{$gene}}[4] = $len;
			${$hash{$gene}}[5] = $evalue;
		}else{
			if((${$hash{$gene}}[4] < $len) or (${$hash{$gene}}[4]==$len and ${$hash{$gene}}[3]<$identity)){
				${$hash{$gene}}[3] = $identity;
				${$hash{$gene}}[4] = $len;
				${$hash{$gene}}[5] = $evalue;
			}
		}
	}else{
		$hash{$gene} = ["NA", "NA", "NA", $identity, $len, $evalue];
	}
}

print TGT "Gene_ID\tOldRun_identity\tOldRun_len\tOldRun_evalue\tNewRun_identity\tNewRun_len\tNewRun_evalue\n";
foreach my $key (sort(keys %hash)){
	print TGT $key, "\t", join("\t", @{$hash{$key}}), "\n";
}

close OLD;
close NEW;
close TGT;
