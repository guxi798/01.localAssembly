#!/usr/bin/perl -w

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $pattern = shift @ARGV;
my $times = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my %hash = ();
while(my $line = <SRC>){
	if($line =~ /^@/){
		my $preline = $line;
		$line = <SRC>;
		chomp $line;
		if($line =~ /$pattern{$times,}/){
			<SRC>;
			<SRC>;
			next;
		}
		print TGT "$preline";
		print TGT "$line\n";
		$line = <SRC>; print TGT "$line";
		$line = <SRC>; print TGT "$line";
	}
}

close SRC;
close TGT;
