#!/bin/perl -w

#use strict;

if (@ARGV<3) {
	die "$0:Not enough arguments\n";
}

my $frblat=$ARGV[0];
#my $fwmatch=$ARGV[1];
my $frseq=$ARGV[1];
my $fwseq=$ARGV[2];
my $line;
my $seq;
my @temp;
my @rna1;
my @rna2;
my $id;
my $name1;
my $name2;

my $count1=0;
my $count2=0;
my $mark=0;
my $i=0;
my $test=0;

#test the global arguments.
my $len=@ARGV;
print $len,"\n";
#print $fwseq,"\n";
print "Blat file:$ARGV[0];\nFastq file=$ARGV[1];\nCleaned file=$ARGV[2]\n";

#The first round, to pick up the pvalue >=90.

#Build a hash to keep the sequence ID

my %hash=();

open(FRB,"$frblat")||die "Can not open the blat file: $!.\n";
#open(FWM,">$fwmatch")||die "Can not create the highly matched file: $!.\n";
	
foreach $line (<FRB>) {
	chomp $line;
	@temp=split(/\t/,$line);
	$name1=$temp[0];
	my @seq_id=split('_',$name1);	
	my $id=$seq_id[0];
	my $len=$seq_id[1];
	my $len_hit=$temp[2]*$temp[3]/100;

    if ($len_hit/$len>=0.9) {
	    $hash{$id}=1;
    }#first if
}#first foreach

#print FWM "The total number of matches are: $count1.\n";


close(FRB);
#close(FWM);

#The second round, to remove the matched sequence.
$count1=0;
#open(FRM,"$fwmatch")||die "Can not open the highly matched file: $!.\n";
open(FRS,"$frseq")||die "Can not open the RNAseq file: $!.\n";
open(FWS,">$fwseq")||die "Can not create the cleaned RNAseq file: $!.\n";


my $fq_num=0;
my $fq_out_num=0;
my $mark1 = 0;
while (my $line = <FRS>) {
	if($line=~ m/^>/) {
		$fq_num++;
		@temp=split(/\s+/,$seq);
		my $id = $temp[0];
		$id =~ s/^>//;
		if(not exists $hash{$id}){
			$mark1 = 1;
			$fq_out_num++;
			print FWS "$line";
		}else{
			$mark1 = 0;
		}
	}#if
	
	if($mark1){print FWS "$line";}
}

my @key=keys %hash;
my $ncrna_num=@key+0;

print "total reads: $fq_num\n ncRNA reads: $ncrna_num\n left reads: $fq_out_num\n";
        
        
			 
			 