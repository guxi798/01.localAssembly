#!/usr/bin/perl -w
# run the script: time perl 00.script/c3.combine.fasta.pl 06.assembly/03.bowtie.nucl/run.0

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;

## start running the script
opendir(SRC, $srcfolder) or die "Cannot open $srcfolder: $!";
my @subs = sort(grep(/[0-9]+/, readdir(SRC)));
my $file = "Trinity.all.fasta";
system(":> $srcfolder/$file");

foreach my $sub (@subs){
	system("cat $srcfolder/$sub/Trinity.fasta >> $srcfolder/$file");
}

close SRC;
