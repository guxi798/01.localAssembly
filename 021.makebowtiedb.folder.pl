#!/usr/bin/perl -w
# run the script: time perl 00.script/021.makebowtiedb.folder.pl 01.data/05.splitGenes/02.Transcript/run.1 Sapelo 10

use strict;
system("echo 'Running 021.makebowtieb.folder.pl ....' >> job.monitor.txt");

my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $thread = 1;

## check if previous step has succesfully finished
while(1){
	if(!(-e "00.script/shell.script/transfer.saturate.seq.log")){
		sleep $sleeptime; # job hasn't finished, no log file
	}else{ # job has finished
		if(!-s "00.script/shell.script/summary.error.log"){ # check if all jobs are successful
			print "All jobs have been finished successfully\n"; # error file is empty
			last;
		}else{
			die "ERROR: something went wrong in previous steps\n";	
		}
	}
}

## start running the script
opendir(SRC, $tgtfolder) or die "ERROR: Cannot open $tgtfolder: $!";
my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));

system("rm -rf 00.script/shell.script.previous");
system("mv 00.script/shell.script 00.script/shell.script.previous");
system("mkdir -p 00.script/shell.script");

foreach my $sub (@subs){
	my $shell = "00.script/shell.script/makebowtiedb.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N 00.script/shell.script/makebowtiedb.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=1:00:00\n";
	    print SHL "#PBS -l mem=2gb\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

    print SHL "\n";
    print SHL "cd \$PBS_O_WORKDIR\n";
	if($platform eq "sapelo"){
    	print SHL "module load bowtie2/2.2.4\n";
	}elsif($platform eq "zcluster"){
		print SHL "PATH=/usr/local/bowtie2/2.2.3/bin/:\$PATH";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	print SHL "cd $tgtfolder/$sub\n";
	print SHL "time bowtie2-build -f -q $sub.fasta $sub\n";
	print SHL "cd ../../../../../\n";
	
	close(SHL);
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "zcluster"){
		system("qsub -q rcc-30d -pe thread $thread $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

close(SRC);
system("echo 'Finished 021.makeblastdb.folder.pl!' >> job.monitor.txt");
