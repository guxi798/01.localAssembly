#!/usr/bin/perl -w
# run the script: time perl 00.script/06.assembly.trinity.folder.pl 04.retrieve.reads/03.bowtie.nucl/run.0 06.assembly/03.bowtie.nucl/run.0 genome Sapelo 1

use strict;
system("echo 'Running 06.assembly.trinity.folder.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $tgtfolder =~ /run\.([0-9]+)/;
my $thread = 8;
my $memory = 80;

## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));

=pod
if($run == 0 and $mode eq 'genome'){
	my @temp = @chks;
	my @R1 = ();
	my @R2 = ();
	my @R3 = ();
	foreach my $i (@temp){
		push @R1, "$srcfolder/$i/unmap.$i.R1.fasta";
		push @R2, "$srcfolder/$i/unmap.$i.R2.fasta";
		if(-e "$srcfolder/$i/unmap.$i.long.fasta"){
			push @R3, "$srcfolder/$i/unmap.$i.long.fasta";
		}
	}
	
	push @chks, 99999;
	system("mkdir -p $srcfolder/99999");
	system("cat ".join(" ", @R1)." > $srcfolder/99999/retrieved.99999.R1.fasta");
	system("cat ".join(" ", @R2)." > $srcfolder/99999/retrieved.99999.R2.fasta");
	if(scalar @R3 > 0){
		system("cat ".join(" ", @R3)." > $srcfolder/99999/retrieved.99999.long.fasta");
	}
}
=cut
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("retrievebowtie.reads.$chk.o*");
		my @stdout = glob("retrievebowtie.reads.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			#print "$stderr[0]\n$stdout[0]\n";
			system("grep -E 'ERROR|Error|error' retrievebowtie.reads.$chk.e* >> 00.script/04.retrieve.script/run.$run/summary.error.log\n");
			#system("echo 'success' > 00.script/shell.script/retrievebowtie.reads.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/04.retrieve.script/run.$run/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/retrieved.$chk.R1.fasta" and -s "$srcfolder/$chk/retrieved.$chk.R2.fasta")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					system("rm -f retrievebowtie.reads.$chk.*");
					system("echo 'Resubmitting the job: retrievebowtie.reads.$chk.sh' >> job.monitor.txt");
					system("qsub 00.script/04.retrieve.script/run.$run/retrievebowtie.reads.$chk.sh");
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

system("mv retrievebowtie.reads.* 00.script/04.retrieve.script/run.$run/");
system("chmod 777 -R 00.script/04.retrieve.script/run.$run");
system("chmod 777 -R $srcfolder");

system("rm -rf 00.script/06.trinity.script/run.$run");
system("mkdir -p 00.script/06.trinity.script/run.$run");

foreach my $sub (@subs){
	my $shell = "00.script/06.trinity.script/run.$run/assembly.trinity.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N assembly.trinity.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=120:00:00\n";
	    print SHL "#PBS -l mem=",$memory,"gb\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

    print SHL "\n";
    print SHL "cd \$PBS_O_WORKDIR\n";
    my $command1 = 0;
	if($platform eq "sapelo"){
		$command1 = "/usr/local/apps/trinity/r20140717";
	    print SHL "module load trinity/r20140717\n\n";
    }elsif($platform eq "zcluster"){
		$command1 = "/usr/local/trinity/r20140717";
		print SHL "export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:\${LD_LIBRARY_PATH}\n";
		print SHL "export PATH=/usr/local/gmap-gsnap/latest/bin/:\${PATH}\n\n";
		$thread /= 2;
		$memory = 25;
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	if(-e "$srcfolder/$sub/retrieved.$sub.long.fasta"){
		print SHL "time $command1/Trinity --seqType fa --CPU $thread --JM ",$memory,"G --left $srcfolder/$sub/retrieved.$sub.R1.fasta --right $srcfolder/$sub/retrieved.$sub.R2.fasta --long_reads $srcfolder/$sub/retrieved.$sub.long.fasta --output $tgtfolder/$sub --min_contig_length 100 \n";
	}else{
		print SHL "time $command1/Trinity --seqType fa --CPU $thread --JM ",$memory,"G --left $srcfolder/$sub/retrieved.$sub.R1.fasta --right $srcfolder/$sub/retrieved.$sub.R2.fasta --output $tgtfolder/$sub --min_contig_length 100 \n";
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
system("echo 'Finished 06.assembly.trinity.folder.pl!' >> job.monitor.txt");

system("chmod 777 -R 00.script/06.trinity.script/run.$run");
