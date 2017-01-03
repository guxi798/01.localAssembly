#!/usr/bin/perl -w
# run the script: time perl 00.script/10.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.0 07.map.back/03.bowtie.nucl/run.0 01.data/05.splitGenes/02.Transcript/run.1 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 

use strict;
use Bio::SeqIO;

system("echo 'Running 10.transfer.saturate.seq.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $blastfolder = shift @ARGV;			## blast result folder
my $srcfolder = shift @ARGV;			## assembled contig folder
my $dbfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;			## target folder to put contigs for new run
my $mode = shift @ARGV;					## abs: absolute value; pct: percent value
my $cutoff = shift @ARGV;				## absolute AA number, or percent

## record file: error and output file
system("mkdir -p $tgtfolder");
system(":> $tgtfolder/full.length.contigs.nucl.fasta");
system(":> $tgtfolder/full.length.contigs.prot.fasta");

## read in blast result file
opendir(SRC, $blastfolder) or die "ERROR: Cannot open $blastfolder: $!";
my @subs = sort(grep(/^[0-9]+$/, readdir SRC));
closedir SRC;

foreach my $sub(@subs){
	system("mkdir -p $tgtfolder/$sub");
	system("time perl 00.script/10.detect.full.length.seq.pl $blastfolder/$sub/$sub.contigs.blast.out $srcfolder/$sub/Trinity.new.fasta $dbfolder/$sub/$sub.fasta $tgtfolder/$sub $sub $mode $cutoff > $tgtfolder/$sub/detect.full.length.log");
	system("cat $tgtfolder/$sub/$sub.full.length.contigs.nucl.fasta >> $tgtfolder/full.length.contigs.nucl.fasta");
	system("cat $tgtfolder/$sub/$sub.full.length.contigs.prot.fasta >> $tgtfolder/full.length.contigs.prot.fasta");
}

system("grep -E 'ERROR|Error|error' 00.script/shell.script/transfer.saturate.seq.e > 00.script/shell.script/summary.error.log");
system("echo 'success' > 00.script/shell.script/transfer.saturate.seq.log");

system("echo 'Finished 10.transfer.saturate.seq.pl!' >> job.monitor.txt");

