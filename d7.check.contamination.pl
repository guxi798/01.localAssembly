#!/usr/bin/perl -w

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

my %hash = ();
while(my $line = <SRC>){
	if($line =~ /^@/){
		$line = <SRC>;
		chomp $line;
		if(not exists $hash{$line}){
			$hash{$line} = 1;
		}else{
			$hash{$line} ++;
		}
		<SRC>;
		<SRC>;
	}
}

foreach my $name (sort {$hash{$b} <=> $hash{$a}} keys %hash){
	print TGT  "$name\t$hash{$name}\n";
}

close SRC;
close TGT;
