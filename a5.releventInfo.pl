#!/usr/bin/perl -w
# run the script: time perl 00.script/a5.releventInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/00.PriorData/Ptrichocarpa_210_v3.0.gene.gff3 01.data/04.GeneOfInterest/GeneID.v1.txt 1000

use strict;

my $reffile = shift @ARGV;
my $profile = shift @ARGV;
my $trafile = shift @ARGV;
my $genefile = shift @ARGV;
my $newfile = shift @ARGV;
my $splitSize = shift @ARGV;

open(REF, $reffile)  or die "Cannot find $reffile: $!\n";
open(PRO, $profile) or die "Cannot find $profile: $!\n";
open(TRA, $trafile) or die "Cannot find $trafile: $!\n";
open(GEN, $genefile) or die "Cannot find $genefile: $!\n";
open(TGT, ">$newfile") or die "Cannot find $newfile: $!\n";

## put GOI into hash structure for reference
my %hash = ();
foreach my $line (<REF>){
	chomp $line;
	if($line eq ""){next;}
	my @lines = split(/\t/, $line);
	if(not exists $hash{$lines[0]}){
		#$hash{$lines[0]} = [$lines[1]];
		$hash{$lines[0]} = ();
	}
}
close REF;

## extract protein length info and put into hash
my $len = 0;
my $pre = 0;
my $gene = 0;
while(my $line = <PRO>){
	chomp $line;
	if($line =~ /^>/){
		$line =~ s/>//;
		$gene = $line;
		if(exists $hash{$pre}){
			push @{$hash{$pre}}, $len; 
		}
		$len = 0;
	}else{
		$len = $len + length($line);
		$pre = $gene;
	}
}
push @{$hash{$pre}}, $len;

## extract transcript length info and put into hash
$len = 0;
$pre = 0;
$gene = 0;
my $gc = 0;
while(my $line = <TRA>){
	chomp $line;
	if($line =~ /^>/){
		$line =~ s/>//;
		$gene = $line;
		if(exists $hash{$pre}){
			push @{$hash{$pre}}, $len; 
			push @{$hash{$pre}}, $gc/$len*100;
		}
		$len = 0;
		$gc = 0;
	}else{
		$len = $len + length($line);
		$gc = $gc + ($line =~ tr/gcGC//);
		$pre = $gene;
	}
}
push @{$hash{$pre}}, $len;
push @{$hash{$pre}}, $gc/$len*100;

## extract gene length info and put into hash
<GEN>;
<GEN>;
while(my $line = <GEN>){
	chomp $line;
	my @lines = split(/\t/, $line);
	next if($lines[2] ne "gene");
	my ($start, $end) = ($lines[3], $lines[4]);
	($gene) = $lines[8] =~ /Name=([\w|\.]+)/;
	if(exists $hash{$gene}){
		push @{$hash{$gene}}, abs($start-$end);
	}
}

## write to output file
print TGT "Gene_ID\tProtein_len\tTranscript_len\tGC_Content\tGene_len\tGroup_ID\n";
open(REF, $reffile);
my $i = 1;
my $group = 0;
foreach my $line (<REF>){
	chomp $line;
	if($line eq ""){
		$i++;
		next;
	}
	my @lines = split(/\t/, $line);
    #print "$lines[0]\t",scalar @{$hash{$lines[0]}}, "\n";
	print TGT $lines[0], "\t", join("\t", @{$hash{$lines[0]}}), "\t";
	
	$group = $i * $splitSize;
	print TGT $group, "\n";
}

close REF;
close PRO;
close TRA;
close GEN;
close TGT;
