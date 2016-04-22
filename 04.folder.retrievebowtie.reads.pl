#!/usr/bin/perl -w
# run the script: time perl 00.script/04.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.0 04.retrieve.reads/03.bowtie.nucl/run.0 nucl bowtie.log/bowtie.run.0 Sapelo 10

use strict;
system("echo 'Running 04.folder.retrievebowtie.reads.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $seqtype = shift @ARGV;
my $logfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $tgtfolder =~ /run\.([0-9]+)/;
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
		my @stderr = glob("bowtie.$seqtype.$chk.o*");
		my @stdout = glob("bowtie.$seqtype.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' bowtie.$seqtype.$chk.e* >> 00.script/$logfolder/summary.error.log\n");
			#system("echo 'success' > 00.script/$logfolder/bowtie.$seqtype.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/$logfolder/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/bowtie.out.$chk.sam")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					system("rm -f bowtie.$seqtype.$chk.*");
					system("echo 'Resubmitting the job: bowtie.$seqtype.$chk.sh' >> job.monitor.txt");
					system("qsub 00.script/$logfolder/bowtie.$seqtype.$chk.sh");
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

system("mv bowtie.* 00.script/$logfolder");
system("chmod 777 -R 00.script/$logfolder");
system("chmod 777 -R $srcfolder");

system("rm -rf 00.script/04.retrieve.script/run.$run");
system("mkdir -p 00.script/04.retrieve.script/run.$run");

foreach my $sub (@subs){
	my $shell = "00.script/04.retrieve.script/run.$run/retrievebowtie.reads.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
 	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N retrievebowtie.reads.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=12:00:00\n";
	    print SHL "#PBS -l mem=40gb\n";
		print SHL "\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	print SHL "mkdir -p $tgtfolder/$sub\n";
	
	my $unmap = "map";
	if($run == 0){
		$unmap = "both";
	}
	print SHL "time perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.sam $unmap $tgtfolder/$sub/retrieved.$sub.R1.fasta $tgtfolder/$sub/unmap.$sub.R1.fasta $tgtfolder/$sub/retrieved.$sub.R2.fasta $tgtfolder/$sub/unmap.$sub.R2.fasta\n";
	if(-e "$srcfolder/$sub/bowtie.out.$sub.long.sam" and -s "$srcfolder/$sub/bowtie.out.$sub.long.sam"){
		print SHL "time perl 00.script/04.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.long.sam $unmap $tgtfolder/$sub/retrieved.$sub.long.fasta $tgtfolder/$sub/unmap.$sub.long.fasta\n";
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

close SRC;
system("echo 'Finished 04.folder.retrievebowtie.reads.pl!' >> job.monitor.txt");

system("chmod 777 -R 00.script/04.retrieve.script/run.$run");
