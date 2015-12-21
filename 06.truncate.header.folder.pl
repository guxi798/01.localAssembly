#!/usr/bin/perl -w
# run the script: time perl 00.script/06.truncate.header.folder.pl 06.assembly/03.bowtie.nucl/run.0 Sapelo 60 

use strict;
system("echo 'Running 06.truncate.header.folder.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $srcfolder =~ /(run\.[0-9]+)/;
my $thread = 1;

## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
#if($run eq 'run.0'){
	#push @chks, 99999;
#}
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	#print "I'm here\n";
	while(my $chk = shift @temp){
		my @stderr = glob("assembly.trinity.$chk.o*");
		my @stdout = glob("assembly.trinity.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' assembly.trinity.$chk.e* >> 00.script/trinity.script/$run/summary.error.log\n");
			#system("echo 'success' > 00.script/shell.script/assembly.trinity.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/trinity.script/$run/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-e "$srcfolder/$chk/Trinity.fasta")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f assembly.trinity.$chk.*");
					#system("rm -r $srcfolder/$chk");
					system("echo 'Resubmitting the job: assembly.trinity.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/trinity.script/$run/assembly.trinity.$chk.sh");
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

system("mv assembly.trinity.* 00.script/trinity.script/$run/");
system("rm -rf 00.script/shell.script.previous");
system("mv 00.script/shell.script 00.script/shell.script.previous");
system("mkdir -p 00.script/shell.script");

foreach my $sub (@subs){
	if($sub eq '99999'){
		next;
	}
	my $shell = "00.script/shell.script/truncate.header.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N 00.script/shell.script/truncate.header.$sub\n";
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

	print SHL "time perl 00.script/06.truncate.header.pl $srcfolder/$sub/Trinity.fasta $srcfolder/$sub/Trinity.new.fasta\n\n";
	
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
system("echo 'Finished 06.truncate.header.folder.pl!' >> job.monitor.txt");

