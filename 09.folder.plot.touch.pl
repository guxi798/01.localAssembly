#!/usr/bin/perl -w
# run the script: time perl 00.script/09.folder.plot.touch.pl 07.map.back/03.bowtie.nucl/run.0 09.plot/03.bowtie.nucl/run.0 Sapelo 10

use strict;
system("echo 'Running 09.folder.plot.touch.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $thread = 1;

## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("blastx.back.$chk.o*");
		my @stdout = glob("blastx.back.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' blastx.back.$chk.e* >> 00.script/shell.script/summary.error.log\n");
			#system("echo 'success' > 00.script/shell.script/summarize.blastx.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/shell.script/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/$chk.contigs.blast.out")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					system("rm -f blastx.back.$chk.*");
					system("echo 'Resubmitting the job: blastx.back.$chk.sh' >> job.monitor.txt");
					system("qsub 00.script/shell.script/blastx.back.$chk.sh");
					#last; # There are some jobs failed
				}
				else{
					$i++;
				}
			}
			if($i == scalar @chks){
				system("echo 'All jobs have been finished successfully' >> job.monitor.txt"); # error file is empty
				last;
			}
		}else{
			die "ERROR: something went wrong in previous steps\n";	
		}
	}
	sleep $sleeptime;
}
close CHK;

## start running the script
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));

system("mv blastx.back.* 00.script/shell.script/");
system("rm -rf 00.script/shell.script.previous");
system("mv 00.script/shell.script 00.script/shell.script.previous");
system("mkdir -p 00.script/shell.script");

foreach my $sub (@subs){
	my $shell = "00.script/shell.script/plot.blastx.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N plot.blastx.$sub\n";
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

    print SHL "mkdir -p $tgtfolder/$sub\n";
    print SHL "time echo 'test' > $tgtfolder/$sub/$sub.evalue.emf\n";
    print SHL "time echo 'test' > $tgtfolder/$sub/$sub.identity.emf\n";
    print SHL "time echo 'test' > $tgtfolder/$sub/$sub.length.emf\n";
	
	print SHL "grep -E 'ERROR|Error|error' plot.blastx.$sub.e >> 00.script/shell.script/summary.error.log\n";
	print SHL "echo 'success' > 00.script/shell.script/plot.blastx.$sub.log\n";
	
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
system("echo 'Finished 09.folder.plot.touch.pl!' >> job.monitor.txt");
