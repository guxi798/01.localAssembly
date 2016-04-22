#!/usr/bin/perl -w
# run the script: time perl ../../00.script/d15.folder.get.unplant.seq.pl blastn.batch transcriptome.part3.v3.fa 0.8 zcluster

use strict;
use Bio::SeqIO;
use Bio::SearchIO;

my $reffolder = shift @ARGV;
my $srcfile = shift @ARGV;
my $cut = shift @ARGV;
my $platform = lc(shift @ARGV);

opendir(REF, $reffolder);
my @files = sort(grep(/^rh/, readdir(REF)));
closedir REF;

foreach my $file (@files){
	print "Processing $file...\n";
	my $shell = "$file.clear.sh";
	open(SHL, ">$shell");
	
	if($platform eq 'sapelo'){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N $file.clear.\n";
		print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
		print SHL "#PBS -l walltime=480:00:00\n";
		print SHL "#PBS -l mem=10gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n\n";
		
		print SHL "time ../../00.script/d15.get.unplant.seq.pl $reffolder/$file/$file.xml.out $reffolder/$file/$file.hash.out $cut\n";
		
		close SHL;
		system("chmod u+x $shell");
		system("qsub $shell");
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n\n";
		
		print SHL "time ../../00.script/d15.get.unplant.seq.pl $reffolder/$file/$file.xml.out $reffolder/$file/$file.hash.out $cut\n";
		
		system("chmod u+x $shell");
		system("qsub -q rcc-30d $shell");
	}else{
		close SHL;
		die "Please provide the platform: Sapelo or Zcluster.\n";
	}
}

