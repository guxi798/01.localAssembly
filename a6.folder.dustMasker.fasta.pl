#!/usr/bin/perl -w
# run the script: time perl 00.script/a6.folder.dustMasker.fasta.pl 01.data/02.Fasta
use strict;

my $srcfolder = shift @ARGV;
opendir(SRC, $srcfolder);
my @sams = sort(grep(/^[A-Z][a-z]\w+G\w*/, readdir SRC));
closedir SRC;
my @R = ("R1", "R2");

system("mkdir -p 00.script/dustMask.script");
system("rm -f 00.script/dustMask.script/*");

foreach my $sam (@sams){
	foreach my $R (@R){
		my $shell = "00.script/dustMask.script/dustMask.$sam.$R.sh";
		open(SHL, ">$shell");
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N dustMask.$sam.$R\n";
		print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=20gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n\n";
		print SHL "time /usr/local/apps/ncbiblast+/2.2.29/bin/dustmasker -in $srcfolder/$sam/$sam.$R.fasta -out $srcfolder/$sam/$sam.$R.masked.fasta -outfmt fasta\n";
		close SHL;
		
		system("qsub $shell");
	}
}