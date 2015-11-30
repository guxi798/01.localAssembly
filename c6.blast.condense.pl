#!/usr/bin/perl -w
# run the script: time perl 00.script/c3.combine.fasta.pl 06.assembly/03.bowtie.nucl/run.0
# this condense blast results by combining contigs with non-overlap multiple alignment

use strict;

## read in parameters required by the script
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

## start running the script
open(SRC, $srcfile) or die "Cannot open $srcfile: $!";
open(TGT, ">$tgtfile");

my $i = 0;
my $gene = 0;
my $pregene = 0;
my $start = 0;
my $prestart = 0;
my $end = 0;
my $preend = 0;
my $len = 0;
my $prelen = 0;
my @prelines = ();

my %hash = ();

foreach my $line (<SRC>){
	$i++;
	chomp $line;
	my @lines = split(/\t/, $line);
	$gene = $lines[0];
	$start = $lines[6] < $lines[7]? $lines[6]:$lines[7];
	$end = $lines[6] < $lines[7]? $lines[7]:$lines[6];
	$len = $end - $start + 1;
	$lines[6] = $start;
	$lines[7] = $end;
	$lines[3] = $len;
	
	if(not exists $hash{$gene}){
		$hash{$gene} = [\@lines];
	}else{
		my @arrays = @{$hash{$gene}};
		for(my $i = 0; $i < scalar @arrays; $i++){
			my $array = $arrays[$i];
			@prelines = @$array;
			if($lines[7] < $prelines[6]){
				unshift(@{$hash{$gene}}, \@lines);
				last;
			}elsif(($lines[6] < $prelines[6]) and ($prelines[6] < $lines[7]) and ($lines[7] < $prelines[7])){
				$prelines[6] = $start;
				$prelines[3] = $prelines[7] - $start + 1;
				$@{$hash{$gene}}[$i] = \@prelines;
				last;
			}elsif(($prelines[6] < $lines[6]) and ($lines[6] < $prelines[7]) and ($prelines[7] < $lines[7])){
				$prelines[7] = $end;
				$prelines[3] = $end - $prelines[6] + 1;
				$@{$hash{$gene}}[$i] = \@prelines;
				last;
			}elsif($lines[6] <= $prelines[6] and $prelines[7] <= $lines[7]){
				$@{$hash{$gene}}[$i] = \@lines;
				last;
			}elsif($prelines[6] <= $lines[6] and $lines[7] <= $prelines[7]){
				last;
			}elsif($lines[6] > $prelines[7]){
				if($i == scalar @arrays -1){
					push @{$hash{$gene}}, \@lines;
				}else{
					next;
				}
			}
		}
	}
	
}

foreach my $key (sort keys(%hash)){
	my @arrays = @{$hash{$key}};
	my $len = 0;
	my @lines = ();
	foreach my $array (@arrays){
		@lines =  @$array;
		$len += $lines[3];
	}
	$lines[3] = $len;
	print TGT join("\t", @lines), "\n";
}

close SRC;
close TGT;
