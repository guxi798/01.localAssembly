#!/usr/bin/perl -w
# run the script: time perl ../../00.script/d15.get.unplant.seq.pl blastn.batch/rh1/rh1.xml.out blastn.batch/rh1/rh1.hash.out 0.8

use strict;
use Bio::SeqIO;
use Bio::SearchIO;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $cut = shift @ARGV;

open(TGT, ">$tgtfile");
my %hash = ();
my $blast = new Bio::SearchIO(-format => 'blastxml', -file => "$srcfile");
while(my $result = $blast->next_result){
	my $qlen = $result->query_length;
	while(my $hit = $result->next_hit){
		my $hitname = $hit->name;
		my $hitlen = $hit->length;
		my $len = &min($qlen, $hitlen);
		#my $len = $hitlen;
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		my $frac = $identity/$len;
				if($frac > $cut){
			if(not exists $hash{$hitname}){
				$hash{$hitname} = 0;
				#print "$hitname\t$qlen\t$hitlen\t$identity\t$frac\n";
				print TGT "$hitname\t$qlen\t$hitlen\t$identity\t$frac\n";
			}
		}
	}
}

close TGT;

## function
sub min{
	my $l = shift @_;
	my $r = shift @_;
	if($l <= $r*0.7){return $r;}
	if($l <= $r){return $l}
	else{return $r;}
}