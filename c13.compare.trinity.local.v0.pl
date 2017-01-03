#!/usr/bin/perl -w
# run the script: time 00.script/c13.compare.Trinity.local.pl 05.Trinity/count1 08.full.length/Final.fasta 10.unmapped.reads.trinity/full.length.contigs.fasta

use strict;
use Bio::SeqIO;

my $trinityfile = shift @ARGV;
my $localfile = shift @ARGV;
my $newfile = shift @ARGV;

my $seqio_obj = Bio::SeqIO->new(-file => $localfile, -format => "fasta");
my $count1 = 0;
my $count2 = 0;
my $count3 = 0;
my %hash = ();
while(my $seq_obj = $seqio_obj->next_seq){
	$count1 ++;
	my $line = $seq_obj->id;
	my $desc = $seq_obj->desc;
	my @lines = split(/\s+/, $desc);
	if(not exists $hash{$lines[1]}){
		$hash{$lines[1]} = 0;
		$count2 ++;
	}
}

print "Local: $count1 contigs\n";
print "Local: $count2 unique genes\n";

open(SRC, $trinityfile);
$count1 = 0;
$count2 = 0;
$count3 = 0;
my %trihash = ();
my %triuniq = ();
while(my $line = <SRC>){
	chomp $line;
	$count1 ++;
	my @lines = split(/\s+/, $line);
	if(exists $hash{$lines[1]} and not exists $trihash{$lines[1]}){
		$count2 ++;
		$trihash{$lines[1]} = 0;
	}elsif(not exists $trihash{$lines[1]}){
		$trihash{$lines[1]} = 0;
		$triuniq{$lines[1]} = 0;
		$count3 ++;
	}
}
close SRC;

print "Trinity: $count1 contigs\n";
print "Trinity: ", $count2+$count3, " unique genes\n";
print "Trinity: $count3 not shared by local assembly\n";

$seqio_obj = Bio::SeqIO->new(-file => $newfile, -format => "fasta");
$count1 = 0;
$count2 = 0;
my %newhash = ();
my %hash1 = ();
while(my $seq_obj = $seqio_obj->next_seq){
	$count1 ++;
	my $line = $seq_obj->id;
	my $desc = $seq_obj->desc;
	my @lines = split(/\s+/, $desc);
	if(exists $triuniq{$lines[1]} and not exists $hash1{$lines[1]}){
		$hash1{$lines[1]} = 0;
		$count2 ++;
	}elsif(not exists $newhash{$lines[1]} and not exists $hash{$lines[1]}){
		$newhash{$lines[1]} = 0;
		$count3 ++;
	}
}

print "New recover: $count1 contigs\n";
print "New recover: share $count2 Trinity unique genes\n";
print "New recover: give $count3 more genes\n";