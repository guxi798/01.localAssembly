#!/usr/bin/perl -w
# run the script: time perl 00.script/c14.folder.bowtie.pl 01.data/00.PriorData/StHeMenE 01.data/01.Fastq 14.bowtie.aligment paired-end Sapelo

use strict;

## read in parameters required by the script
my $reffile = shift @ARGV;				## reference file with gene ID and protein length info
my $srcfolder = shift @ARGV;			## reads file
my $tgtfolder = shift @ARGV;			## target folder to put contigs for new run
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $thread = 1;
my $scriptfolder = "00.script/14.bowtie.script";

## read in blast result file
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^\w+/, readdir SRC));
closedir SRC;

system("rm -rf $scriptfolder");
system("mkdir -p $scriptfolder");

foreach my $sub (@subs){	## loop over parallized groups
	my $file = 0;
	if($sub =~ /F$|Fu$|R$/){
		$mode = 'single-end'; 
		$file = "$sub.fna";
	}
 	my $shell = "$scriptfolder/bowtie.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N bowtie.$sub\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=15gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	my $t = $thread / 2;
		if($platform eq "sapelo"){
    	print SHL "module load bowtie2/2.2.4\n";
	}elsif($platform eq "zcluster"){
		print SHL "PATH=/usr/local/bowtie2/2.2.3/bin/:\$PATH\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	if($mode eq "paired-end"){
		print SHL "mkdir -p $tgtfolder/$sub\n";
		print SHL "time bowtie2 -x $reffile -p $thread -1 $srcfolder/$sub/$sub.R1.fastq -2 $srcfolder/$sub/$sub.R2.fastq -S $tgtfolder/$sub/$sub.sam\n";
	}elsif($mode eq "single-end"){
		print SHL "mkdir -p $tgtfolder/$sub\n";
		if($sub =~ /F$|Fu$|R$/){
			print SHL "time bowtie2 -f -x $reffile -p $thread -U $srcfolder/$sub/$file -S $tgtfolder/$sub/$sub.sam\n";
		}else{
			print SHL "time bowtie2 -x $reffile -p $thread -U $srcfolder/$sub/$sub.fastq -S $tgtfolder/$sub/$sub.sam\n";
		}
	}
	
	close SHL;
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "zcluster"){
		system("qsub -q rcc-30d -pe thread $thread $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

}

