#!/usr/bin/perl -w
# after running 10.transfer.sature.seq.pl
# run the script: time perl 00.script/c12.filter.by.length.pl 01.data/06.TargetTranscriptome/transcriptome.part3.v1.fa 01.data/06.TargetTranscriptome/transcriptome.part3.v2.fa 200

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $cut = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");
my $mark = 0;
my $count = 0;
my $count1 = 0;
my $count2 = 0;

while(my $line= <SRC>){
	if($line =~ /^>/){
		$count ++;
		$line =~ s/^>//;
		my @lines = split(/\s+/, $line);
		my $len = $lines[1];
		$len =~ s/len=//;
		if($len >= $cut){
			$count1 ++;
			$mark = 1;
			print TGT ">$line";
		}else{
			$count2 ++;
			$mark = 0;
		}
	}else{
		if($mark){print TGT "$line"}
	}
}
close SRC;
close TGT;

print "There are $count contigs in total\n";
print "There are $count1 contigs with length larger than $cut bp\n";
print "There are $count2 contigs with length shorter than $cut bp\n";