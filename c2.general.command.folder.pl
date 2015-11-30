#!/usr/bin/perl -w
# run the script: time perl 00.script/c2.general.command.folder.pl 01.data/01.Fastq

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;

## start running the script
opendir(SRC, $srcfolder) or die "Cannot open $srcfolder: $!";
my @subs = sort(grep(/\w+/, readdir(SRC)));

foreach my $sub (@subs){
	opendir(SUB, "$srcfolder/$sub");
	my @new = sort(grep(/[0-9]+/, readdir(SUB)));
	foreach my $new (@new){
		system("rm -f $srcfolder/$sub/$new/*.fa");
		system("rm -f $srcfolder/$sub/$new/*_count");
		system("rm -f $srcfolder/$sub/$new/*.bam");
		system("rm -f $srcfolder/$sub/$new/*.finished");
		system("rm -f $srcfolder/$sub/$new/*.sam");
		system("rm -f $srcfolder/$sub/$new/*.txt");
		system("rm -f $srcfolder/$sub/$new/*.histo");
		system("rm -f $srcfolder/$sub/$new/*.ebwt");
		system("rm -f $srcfolder/$sub/$new/*.thread_*");
		system("rm -rf $srcfolder/$sub/$new/chrysalis");
	}
	close SUB;
}

close SRC;
