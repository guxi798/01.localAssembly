#!/usr/bin/perl -w

use strict;

my $input1 = shift @ARGV;
my $input2 = shift @ARGV;
my $ref = shift @ARGV;

open(INT1, "$input1");
open(INT2, "$input2");
open(REF, "$ref");

my %hash1 = ();
my %hash2 = ();
foreach my $line(<INT1>){
	chomp $line;
	if(not exists $hash1{$line}){
		$hash1{$line} = 0;
	}
}
foreach my $line(<INT2>){
	chomp $line;
	if(not exists $hash2{$line} and not exists $hash1{$line}){
		$hash2{$line} = 0;
	}
}

my $pre = 0;
my ($sim1, $sim2) = 0;
my ($count1, $count2) = 0;
foreach my $line (<REF>){
	chomp $line;
	my @lines = split(/\t/, $line);
	my $gene = $lines[1];
	my @subs = split(/\./, $gene);
	pop @subs;
	$gene = join(".", @subs);
	#if($gene eq $pre){next;}
	if(exists $hash1{$gene}){
		$count1 ++;
		$sim1 += $lines[2];
	}elsif(exists $hash2{$gene}){
		$count2 ++;
		$sim2 += $lines[2];
	}
}

print "Unique gene of Trinity, average similarity: ", $sim1/$count1, "\n";
print "Common gene of Trinity, average similarity: ", $sim2/$count2, "\n";

close INT1;
close INT2;
close REF;