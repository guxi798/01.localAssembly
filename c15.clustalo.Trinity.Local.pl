#!/usr/bin/perl -w
#run the script: time perl 00.script/c15.clustalo.Trinity.Local.pl 13.Local.Trinity.alignment/01.data 13.Local.Trinity.alignment/02.alignment

use strict;

my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;

system("mkdir -p $tgtfolder");

opendir(SRC, $srcfolder);
my @files = grep(/[0-9]+/, readdir(SRC));
closedir(SRC);

my %groups = ();
foreach my $file (@files){
	my ($temp) = $file =~ /^([0-9]+)/;
	if(not exists $groups{$temp}){
		$groups{$temp} = 0;
	}
}

foreach my $group (sort {$a <=> $b} keys %groups){
	print "$group....\n";
	system("/usr/local/clustal-omega/1.2.1/bin/clustalo -t RNA -i $srcfolder/$group.nucl.fasta -o $tgtfolder/$group.nucl.aln --outfmt=clu --force");
	system("/usr/local/clustal-omega/1.2.1/bin/clustalo -t Protein -i $srcfolder/$group.prot.fasta -o $tgtfolder/$group.prot.aln --outfmt=clu --force");
}