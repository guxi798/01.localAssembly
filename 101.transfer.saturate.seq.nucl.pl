#!/usr/bin/perl -w
# run the script: time perl 00.script/101.transfer.saturate.seq.pl 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/blastx.out 10.unmapped.reads.trinity 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 

use strict;

use Bio::SeqIO;
use Bio::SearchIO;

## read in parameters required by the script
my $usage = "usage: $0 blast+.outfmt6 query.fasta db.fasta output_folder [output_prefix=NameOfBlastFileHere] [mode=abs|pct] [cutoff=NumericCutoff]\n\n";

my $blastfile = shift @ARGV or die $usage;              ## blast result
my $srcfile = shift @ARGV or die $usage;				## assembled contig file
my $dbfile = shift @ARGV or die $usage;
my $tgtfolder = shift @ARGV or die $usage;              ## target folder
my $output_prefix = shift @ARGV || "$blastfile";
my $mode = shift @ARGV || "pct";					## abs: absolute value; pct: percent value
my $cutoff = shift @ARGV || "0.02";                 ## absolute AA number, or percent

system("mkdir -p $tgtfolder");
my $logfile = "$tgtfolder/detect.full.length.log";

open(TGT1, ">$tgtfolder/$output_prefix.full.contigs.nucl.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.full.length.contigs.fasta: $!";
open(TGT2, ">$tgtfolder/$output_prefix.incomplete.contigs.nucl.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.incomplete.fasta: $!";
open(TGT3, ">$tgtfolder/$output_prefix.unmapped.contigs.nucl.fasta");

## read in blast results
my %full_length = ();
my %incomplete_length = ();

my $blast = new Bio::SearchIO(-format 	=> 'blastxml', -file	=> $blastfile);
while(my $result = $blast->next_result){
	my @names = split(/\s+/, $result->query_description);
	my $qname = $names[0];
	my $qlen = $result->query_length;
	while(my $hit = $result->next_hit){
		my $hitname = $hit->name;
		my $hitlen = $hit->length;
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		my $frac_iden1 = $hit->frac_identical;
		my $len = 0;
		while(my $hsp = $hit->next_hsp){
			$len += $hsp->length('hit');
		}
		my $frac_iden2 = $len/$hitlen;
		
		if(($mode eq "abs" and $qlen > ($hitlen-$cutoff)) || ($mode eq "pct" and $frac_iden2 >= 1-$cutoff)){
			if(not exists $full_length{$qname}){
				$full_length{$qname} = { query_id => $qname,
										   hit_id => $hitname,
										   query_len => $qlen,
										   hit_len => $hitlen
										   };
			}
		}else{
			if(not exists $incomplete_length{$qname}){
				$incomplete_length{$qname} = { query_id => $qname,
										   hit_id => $hitname,
										   query_len => $qlen,
										   hit_len => $hitlen
										   };
			}
		}
	}
}

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");

while(my $seq_obj = $seqio_obj->next_seq){
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my @temp = split(/\s+/, $seq_obj->desc);
	my $seqlen = $temp[0];
	
	if(exists $full_length{$id}){
		print TGT1 ">$id len=", $full_length{$id}->{query_len}, " ", 
					$full_length{$id}->{hit_id}, " ",
					$full_length{$id}->{hit_len}, " ", 
					$full_length{$id}->{query_len}, " blast_full 0 \n";
		print TGT1 "$seq\n";
	}elsif(exists $incomplete_length{$id}){
		print TGT2 ">$id len=", $incomplete_length{$id}->{query_len}, " ", 
				$incomplete_length{$id}->{hit_id}, " ",
				$incomplete_length{$id}->{hit_len}, " ", 
				$incomplete_length{$id}->{query_len}, " 0 \n";
		print TGT2 "$seq\n";
	}else{
		print TGT3 ">$id $seqlen\n";
		print TGT3 "$seq\n";
	}
}

close TGT1;
close TGT2;
close TGT3;

## function
sub min{
	my $l = shift @_;
	my $r = shift @_;
	if($l <= $r){return $l}
	else{return $r;}
}
