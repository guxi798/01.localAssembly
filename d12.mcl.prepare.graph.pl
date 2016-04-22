#!/usr/bin/perl -w
# run the script: time /usr/local/mcl/latest/bin/mcl  12.MCL.remove.redundancy/mcl.graph.txt --abc -o 12.MCL.remove.redundancy/mcl.out.txt -I 5.0

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $mode = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my $preq = "";
my $pres = "";
my $score = 0.0;
my $count = 0;
my $average = 0.0;
my $iter = 0;
my $e = 0;

foreach my $line(<SRC>){
	$iter ++;
	chomp $line;
	my @lines = split(/\t/, $line);
	my ($q, $s) = @lines[0..1];
	if($mode eq "ncbi"){
		$e = $lines[10];
	}elsif($mode eq "wu"){
		$e = $lines[2];
	}else{
		die("Available mode: ncbi or wu.");
	}
	
	## take log of evalue
	if($e == 0.0){
		$e = 181;
	}else{
		$e = -&log10($e);
	}
	
	## get gene id from transcript id
	if($q eq $preq and $s eq $pres){
		$score += $e;
		$count += 1;
	}else{
		if($iter == 1){
			$score = $e;
			$count = 1;
			$preq = $q;
			$pres = $s;
			next;
		}
		if($count == 0){
			$average = $score;
		}else{
			$average = $score/$count;
		}
		print TGT  "$preq\t$pres\t$average\n";
		$score = $e;
		$count = 1;
		$preq = $q;
		$pres = $s;
	}
}

if($count == 0){
	$average = $score;
}else{
	$average = $score/$count;
}
print TGT  "$preq\t$pres\t$average\n";

close(SRC);
close(TGT);

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
