#!/usr/bin/perl -w
# run the script: time perl 00.script/d18.annotation.pl 01.data/06.TargetTranscriptome/transcriptome.fa 14.annotation/annotation.txt 14.annotation/Ptr.blastx.out 01.data/00.PriorData/Ptrichocarpa_210_v3.0.annotation_info.txt 14.annotation/Ath.blastx.out 01.data/00.PriorData/Athaliana_167_TAIR10.annotation_info.txt 14.annotation/Uni.blastx.out 01.data/00.PriorData/uniprot_sprot.fasta

use strict;
use Bio::SeqIO;

my $seqfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $blast = shift @ARGV;
my $ref = shift @ARGV;
my $athblast = shift @ARGV;
my $athref = shift @ARGV;
my $uniblast = shift @ARGV;
my $uniref = shift @ARGV;

my %hash = ();

open(BLT, $blast);
while(my $line = <BLT>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}{'Rsp'}){
		$hash{$lines[0]}{'Rsp'} = $lines[1];
	}
}
close BLT;

open(BLT, $athblast);
while(my $line = <BLT>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}{'Ath'}){
		$hash{$lines[0]}{'Ath'} = $lines[1];
	}
}
close BLT;

open(BLT, $uniblast);
while(my $line = <BLT>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}{'Uni'}){
		$hash{$lines[0]}{'Uni'} = $lines[1];
	}
}
close BLT;

my %refhash = ();
open(REF, $ref); <REF>;
while(my $line = <REF>){
	my @lines = split(/\t/, $line);
	chomp $lines[12];
	if(not exists $refhash{$lines[1]}){
		if($lines[9] eq ''){$lines[9] = '.';}else{$lines[9]=~s/GO:GO/GO/g;}
		if($lines[11] eq ''){$lines[11] = '.';}
		if($lines[12] eq ''){$lines[12] = '.';}
		$refhash{$lines[1]} = [$lines[9], $lines[11], $lines[12]];
	}
}
close REF;

open(REF, $athref); <REF>;
while(my $line = <REF>){
	my @lines = split(/\t/, $line);
	chomp $lines[12];
	if(not exists $refhash{$lines[1]}){
		if($lines[9] eq ''){$lines[9] = '.';}else{$lines[9]=~s/GO:GO/GO/g;}
		if($lines[11] eq ''){$lines[11] = '.';}
		if($lines[12] eq ''){$lines[12] = '.';}
		$refhash{$lines[1]} = [$lines[9], $lines[12]];
		#print "$lines[9]\t$lines[11]\t$lines[12]\n";
	}
}
close REF;

my $seqio_obj1 = Bio::SeqIO->new(-file => $uniref, -format => "fasta");
while(my $seq_obj = $seqio_obj1->next_seq){
	my $id = $seq_obj->id;
	my $desc = $seq_obj->desc;
	my ($name) = $desc =~ /(.+) OS=/;
	if(not exists $refhash{$id}){
		$refhash{$id} = $name;
	}
}

my $typelen = 2;
open(TGT, ">$tgtfile");
my $seqio_obj = Bio::SeqIO->new(-file => $seqfile, -format => "fasta");
while(my $seq_obj = $seqio_obj->next_seq){
	my $id = $seq_obj->id;
	if(exists $hash{$id}){
		my $ortho1 = '.';
		my $ortho2 = '.';
		my $ortho3 = '.';
		if(exists $hash{$id}{'Rsp'}){$ortho1 = $hash{$id}{'Rsp'};}
		if(exists $hash{$id}{'Ath'}){$ortho2 = $hash{$id}{'Ath'};}
		if(exists $hash{$id}{'Uni'}){$ortho3 = $hash{$id}{'Uni'};}
		print "$ortho1\t$ortho2\t$ortho3\n";
		my @annot1 = ('.') x ($typelen+1);
		my @annot2 = ('.') x $typelen;
		my $annot3 = '.';
		if(exists $refhash{$ortho1}){@annot1 = @{$refhash{$ortho1}};}
		if(exists $refhash{$ortho2}){@annot2 = @{$refhash{$ortho2}};}
		if(exists $refhash{$ortho3}){$annot3 = $refhash{$ortho3};}
		#if(exists $refhash{$ortho1}){print join("\t", @{$refhash{$ortho1}}),"\n";}
		print TGT "$id\t$ortho1\t",join("\t", @annot1), "\t$ortho2\t", join("\t", @annot2), "\t$annot3\n";
	}else{
		my @annot = ('.') x ($typelen*2 + 1 + 2 + 1);
		print TGT "$id\t", join("\t", @annot), "\n";
	}
}