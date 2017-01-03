#!/usr/bin/perl -w
# after running 10.transfer.sature.seq.pl
# run the script: time perl 00.script/c11.combine.full.length.pl 10.unmapped.reads.trinity/full.length.contigs.blastn.xml.out 10.unmapped.reads.trinity/full.length.contigs.fasta 10.unmapped.reads.trinity/temp.fasta

use strict;
use Bio::SearchIO;
use Bio::SeqIO;
use List::Util qw(sum);

my $blastfile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $mode = shift @ARGV;

my %blastx = ();
if($mode eq 'self'){
	my $blastxfile = shift @ARGV;
	open(BLT, $blastxfile);
	while(my $line = <BLT>){
		chomp $line;
		my @lines = split(/\s+/, $line);
		my $qname = $lines[0];
		my $hitname = $lines[1];
		my $mismatch = $lines[4];
		my $bitscore = $lines[11];
		if(not exists $blastx{$qname} or $blastx{$qname}->{bitscore} < $bitscore){
			$blastx{$qname} = {hitname => $hitname,
								mismatch => $mismatch,
								bitscore => $bitscore};
		}
	}
	close BLT;
}

my $blast = new Bio::SearchIO(-format 	=> 'blastxml',
								-file	=> $blastfile);
	
my %iteration = ();
if($mode eq "self"){
	my $seqio = Bio::SeqIO->new(-file=>$srcfile, -format=>'fasta');
	while(my $seq = $seqio->next_seq){
		my $id = $seq->id;
		my @desc = split(/\s+/, $seq->description);
		if(not exists $iteration{$id}){
			$iteration{$id} = $desc[2];
		}
	}
}

my %hash = (); ## hash table for non-redundant contig ID
my $cut1 = 0.95;
my $cut2 = 0.95;
my $cut3 = 0.9;
my $target = 0;

while(my $result = $blast->next_result){
	my @names = split(/\s+/, $result->query_description);
	my $qname = $names[0];
	my $qlen = $result->query_length;
	while(my $hit = $result->next_hit){
		my $hitname = $hit->name;
		if($qname eq $hitname){next;}
		my $hitlen = $hit->length;
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		my $frac_iden1 = $hit->frac_identical;
		my $len = &min($qlen, $hitlen);
		my $frac_iden2 = $identity/$len;
		my $align_len_query = 0;
		my $align_len_hit = 0;
		my $num_hsp = $hit->num_hsps;
		while(my $hsp = $hit->next_hsp){
			$align_len_query += $hsp->length('query');
			$align_len_hit += $hsp->length('hit');
		}
		my $frac_align = &max($align_len_query/$qlen, $align_len_hit/$hitlen);

		if(($frac_iden1 >= $cut1 and $frac_iden2 >= $cut2) or ($frac_align >= $cut3 and $frac_iden1 >= 0.99)){
			if($mode eq "self"){
				if($iteration{$qname} < $iteration{$hitname}){
					if(exists $blastx{$qname} and exists $blastx{$hitname} and $blastx{$qname}->{hitname} eq $blastx{$hitname}->{hitname}){
						if(($blastx{$hitname}->{bitscore} > $blastx{$qname}->{bitscore}) 
							or ($blastx{$hitname}->{bitscore} == $blastx{$qname}->{bitscore} and $hitlen >= $qlen) 
							or ($blastx{$hitname}->{bitscore} >= $blastx{$qname}->{bitscore}*0.99 and $hitlen*0.95 > $qlen)){
							$target = $qname;
						}else{
							$target = $hitname;
						}
					}else{
						if($frac_iden1 >= 0.99){
							if($blastx{$hitname}->{bitscore} > $blastx{$qname}->{bitscore}*0.99 and $hitlen*0.95 > $qlen){
								$target = $qname;
							}else{
								$target = $hitname;
							}
						}else{
							next;
						}
					}
				}elsif($iteration{$qname} > $iteration{$hitname}){
					if(exists $blastx{$qname} and exists $blastx{$hitname} and $blastx{$qname}->{hitname} eq $blastx{$hitname}->{hitname}){
						if(($blastx{$qname}->{bitscore} > $blastx{$hitname}->{bitscore}) 
							or ($blastx{$qname}->{bitscore} == $blastx{$hitname}->{bitscore} and $qlen >= $hitlen) 
							or ($blastx{$qname}->{bitscore} >= $blastx{$hitname}->{bitscore}*0.99 and $qlen*0.95 > $hitlen)){
							$target = $hitname;
						}else{
							$target = $qname;
						}
					}else{
						if($frac_iden1 >= 0.99){
							if($blastx{$qname}->{bitscore} > $blastx{$hitname}->{bitscore}*0.99 and $qlen*0.95 > $hitlen){
								$target = $hitname;
							}else{
								$target = $qname;
							}
						}else{
							next;
						}					
					}
				}else{
						next;
				}
				
				if(not exists $hash{$target}){
					$hash{$target} = 0;
				}
				print "$qname\t$qlen\t$hitname\t$hitlen\t$target\t$len\t$identity\t$frac_iden1\t$frac_iden2\t$frac_align\n";
			}elsif($mode eq "query"){
				$target = $qname;
				
				if(not exists $hash{$target}){
					$hash{$target} = 0;
				}
				print "$qname\t$qlen\t$hitname\t$hitlen\t$target\t$len\t$identity\t$frac_iden1\t$frac_iden2\t$frac_align\n";
			}
		}
	}
}

my $seqio = Bio::SeqIO->new(-file=>$srcfile, -format=>'fasta');
open(TGT, ">$tgtfile");

my $mark = 0;
my $count1 = 0;
my $count2 = 0;

while(my $seq = $seqio->next_seq){
	$count2++;
	my $contig = $seq->id;
	my $desc = $seq->description;
	if(not exists $hash{$contig}){
		$count1 ++;
		print TGT ">$contig $desc\n";
		print TGT $seq->seq,"\n";
	}
}

print "There are $count2 full contigs from assemblying unmapped reads\n";
print "There are ", $count1, " full contigs not redundant with local assembly\n";

close TGT;

## function
sub min{
	my $l = shift @_;
	my $r = shift @_;
	if($l <= $r){return $l}
	else{return $r;}
}

sub max{
	my $l = shift @_;
	my $r = shift @_;
	if($l <= $r){return $r}
	else{return $l;}
}