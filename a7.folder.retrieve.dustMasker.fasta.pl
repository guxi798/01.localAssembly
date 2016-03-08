#!/usr/bin/perl -w
# run the script: time perl 00.script/a7.folder.retrieve.dustMasker.fasta.pl 01.data/01.Fastq 01.data/02.Fasta 0.5
use strict;

my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $cut = shift @ARGV;
opendir(SRC, $reffolder);
my @sams = sort(grep(/^[A-Z][a-z]\w+G\w*/, readdir SRC));
closedir SRC;
my @R = ("R1", "R2");

system("mkdir -p 00.script/dustMask.script");
system("rm -f 00.script/dustMask.script/*");

foreach my $sam (@sams){
	foreach my $R (@R){
		my $shell = "00.script/dustMask.script/retrieve.dustMask.$sam.$R.sh";
		open(SHL, ">$shell");
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N retrieve.dustMask.$sam.$R\n";
		print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=10gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n\n";
		print SHL "time 00.script/a7.retrieve.dustMasker.fasta.pl $reffolder/$sam/$sam.$R.fastq $srcfolder/$sam/$sam.$R.masked.fasta $reffolder/$sam/$sam.$R.clean.fastq $cut\n";
		close SHL;
		
		system("qsub $shell");
	}
}