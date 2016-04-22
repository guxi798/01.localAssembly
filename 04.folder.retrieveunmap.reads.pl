#!/usr/bin/perl -w
# run the script: time perl 00.script/04.folder.retrieveunmap.reads.pl 01.data/02.Fasta 04.retrieve.reads/03.bowtie.nucl/run.12 0 Sapelo 10

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $long = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $tgtfolder =~ /run\.([0-9]+)/;
my $thread = 1; 

## start running the script

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^[A-Z]+/, readdir(SRC)));

system("rm -rf 00.script/04.retrieve.script/unmap");
system("mkdir -p 00.script/04.retrieve.script/unmap");
system("mkdir -p $tgtfolder/99999");

foreach my $sub (@subs){
	my $shell = "00.script/04.retrieve.script/unmap/retrieveunmap.reads.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
 	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N retrieveunmap.reads.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=12:00:00\n";
	    print SHL "#PBS -l mem=20gb\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

    print SHL "\n";
    print SHL "cd \$PBS_O_WORKDIR\n";
	
	print SHL "time perl 00.script/04.retrieveunmap.reads.pl $sub $run $long\n";
	
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

close SRC;
