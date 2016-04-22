#!/usr/bin/perl -w
# after running 10.transfer.sature.seq.pl
# run the script: time perl 00.script/c11.combine.full.length.pl 10.unmapped.reads.trinity/full.length.contigs.blastn.xml.out 10.unmapped.reads.trinity/full.length.contigs.fasta 10.unmapped.reads.trinity/temp.fasta

use strict;
use strict;
use Bio::SearchIO;
use List::Util qw(sum);

my $blastfile = shift @ARGV;
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

my $blast = new Bio::SearchIO(-format 	=> 'blastxml',
								-file	=> $blastfile);

my %hash = (); ## hash table for non-redundant contig ID
my $cut1 = 0.9;
my $cut2 = 0.8;

while(my $result = $blast->next_result){
	my @names = split(/\s+/, $result->query_description);
	my $qlen = $result->query_length;
	while(my $hit = $result->next_hit){
		my $hitname = $hit->name;
		my $hitlen = $hit->length;
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		my $frac_iden1 = $hit->frac_identical;
		my $len = &min($qlen, $hitlen);
		my $frac_iden2 = $identity/$len;

		if(($frac_iden1 >= $cut1) and ($frac_iden2 >= $cut2)){
			if(not exists $hash{$names[0]}){
				$hash{$names[0]} = 0;
			}
			print "$names[0]\t$hitname\t$len\t$identity\t$frac_iden1\t$frac_iden2\n";
		}
	}
}

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my $mark = 0;
my $count1 = 40000;
my $count2 = 0;
while(my $line = <SRC>){
	if($line =~ /^>/){
		$count2 ++;
		chomp $line;
		my @lines = split(/\s+/, $line);
		$lines[0] =~ s/^>//;
		if(not exists $hash{$lines[0]}){
			$count1 ++;
			#my @ID = split(/_/, $lines[0]);
			#my $ID = "c".$count1."_".$ID[1]."_".$ID[2];
			#print TGT ">$ID $lines[1] 99\n";
			print TGT "$line 99\n";
			$mark = 1;
		}else{
			$mark = 0;
		}
	}else{
		if($mark){print TGT "$line";}
	}
}

print "There are $count2 full contigs from assemblying unmapped reads\n";
print "There are ", $count1 - 40000, " full contigs not redundant with local assembly\n";

close SRC;
close TGT;

## function
sub min{
	my $l = shift @_;
	my $r = shift @_;
	if($l <= $r){return $l}
	else{return $r;}
}