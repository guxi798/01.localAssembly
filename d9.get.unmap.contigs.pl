#!/usr/bin/perl -w
# run the script: 
# module load perl/5.20.2
# time perl 00.script/c7.get.blast.identity.pl 05.Full.Length/both.full.gene.txt 04.blastn/blast.xml.out ../03.Protein.Subset/07.map.back/02.blastn 05.Full.Length/both.full.identity.txt > 05.Full.Length/both.full.identity.out

use strict;
use Bio::SearchIO;
use List::Util qw(sum);

my $srcfile = shift @ARGV;
#my $tgtfile = shift @ARGV;

my $globin = new Bio::SearchIO(-format 	=> 'blastxml',
								-file	=> $srcfile);
my %blast = ();
print "qry_name\thitname\thit_def\tqry_len\thit_len\taligned_len\taligned_frac\n";
while(my $result = $globin->next_result){
	my $qry_name = $result->query_description;
	my $qry_len = $result->query_length;
	while(my $hit = $result->next_hit){
		my $hitname = $hit->name;
		my $hit_def = $hit->hit_description;
		my $hit_len = $hit->length;
		#print "$names[0]\t$len\t$hitname\t";
		#my $identity = scalar($hit->seq_inds('query', 'identical'));
		#print $identity, "\t",$identity/$len,"\n";
		my $map_len = 0;
		while(my $hsp = $hit->next_hsp){
			my $hsp_len = $hsp->length('total');
			if($hsp_len > $map_len){$map_len = $hsp_len;}
		}
		print "$qry_name\t$hitname\t$hit_def\t$qry_len\t$hit_len\t$map_len\t", $map_len/$hit_len, "\n";
	}
}
