#!/usr/bin/perl -w
# run the script: time perl 00.script/a6.split.fasta.pl 01.data/02.Fasta/OrAe1G 100000
use strict;

my $srcfolder = shift @ARGV;
my $size = shift @ARGV;

my @path = split(/\//, $srcfolder);
my $sam = pop @path;
my @R =("R1", "R2");

system("mkdir -p 00.script/repeatMask.script/$sam");
system("rm -f 00.script/repeatMask.script/$sam/*");

foreach my $R (@R){
	my $srcfile = "$sam.$R.fasta";
	open(SRC, "$srcfolder/$srcfile");
	my $j = 1;
	my $tgtfile = "$srcfolder/$sam.$R.part.$j.fasta";
	print "$tgtfile\n";
	open(TGT, ">$tgtfile");

	my $i = 0;
	while(my $line = <SRC>){
		if($line =~ /^>/){
			$i++;
			if($i <= $size){
				print TGT "$line";
				$line = <SRC>;
				print TGT "$line";
			}else{
				# run job on current output
				close TGT;
				my $shell = "00.script/repeatMask.script/$sam/repeatMask.$sam.$R.part.$j.sh";
				open(SHL, ">$shell");
				print SHL "#PBS -S /bin/bash\n";
				print SHL "#PBS -q batch\n";
				print SHL "#PBS -N repeatMask.$sam.$R.part.$j\n";
				print SHL "#PBS -l nodes=1:ppn=2:AMD\n";
				print SHL "#PBS -l walltime=48:00:00\n";
				print SHL "#PBS -l mem=30gb\n\n";
				print SHL "cd \$PBS_O_WORKDIR\n";
				print SHL "module load repeatmasker/4.0.5\n\n";
				print SHL "time /usr/local/apps/repeatmasker/latest/RepeatMasker -noint -e wublast -pa 2 -species arabidopsis -dir $srcfolder $srcfolder/$sam.$R.part.$j.fasta";
				close SHL;
				
				system("qsub $shell");
				
				# start a new partial seq file
				$j ++;
				$tgtfile = "$srcfolder/$sam.$R.part.$j.fasta";
				open(TGT, ">$tgtfile");
				print "$tgtfile\n";
				$i = 1;
				print TGT "$line";
				$line = <SRC>;
				print TGT "$line";
			}
		}
	}
	close TGT;
		my $shell = "00.script/repeatMask.script/$sam/repeatMask.$sam.$R.part.$j.sh";
		open(SHL, ">$shell");
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N repeatMask.$sam.$R.part.$j\n";
		print SHL "#PBS -l nodes=1:ppn=2:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=30gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
		print SHL "module load repeatmasker/4.0.5\n\n";
		print SHL "time /usr/local/apps/repeatmasker/latest/RepeatMasker -noint -e wublast -pa 2 -species arabidopsis -dir $srcfolder $srcfolder/$sam.$R.part.$j.fasta";
		#print SHL "time /usr/local/apps/ncbiblast+/2.2.29/bin/dustmasker -in 01.data/02.Fasta/OrAe1G/OrAe1G.R1.fasta -out OrAe1G.R1/dustMask.test -outfmt fasta\n";
		close SHL;
	
		system("qsub $shell");

	close SRC;
}