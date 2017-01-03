#!/usr/bin/perl -w
#running command: time perl 00.script/03.bowtie.folder.pl 01.data/02.Fasta 01.data/05.splitGenes/02.Transcript nucl paired-end Sapelo 10

use strict;
system("echo 'Running 03.bowtie.folder.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $qryfolder = shift @ARGV;
my $dbfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $seqtype = shift @ARGV;
my $mode = shift @ARGV;
my $logfolder = shift @ARGV;
my $mode1 = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $tgtfolder =~ /run\.([0-9]+)/;
my $thread = 4;

## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("makebowtiedb.$chk.sh.o*");
		my @stdout = glob("makebowtiedb.$chk.sh.e*");
		my @log = glob("00.script/02.makebowtiedb.script/run.$run/$chk.done.*");
		if(!($stderr[0] and $stdout[0] and $log[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' makebowtiedb.$chk.sh.e* >> 00.script/02.makebowtiedb.script/run.$run/summary.error.log\n");
			#system("echo 'success' > 00.script/shell.script/makebowtiedb.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/02.makebowtiedb.script/run.$run/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$dbfolder/$chk/$chk.1.bt2" and -s "$dbfolder/$chk/$chk.rev.1.bt2")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f makebowtiedb.$chk.*");
					#system("rm -f 00.script/02.makebowtiedb.script/run.$run/$chk.done.log");
					system("echo 'Resubmitting the job: makebowtiedb.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/02.makebowtiedb.script/run.$run/makebowtiedb.$chk.sh");
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
opendir(QRY, $qryfolder) or die "ERROR: Cannot open $qryfolder: $!";
my @subs = sort(grep(/^\w.+/, readdir(QRY)));

opendir(DBF, $dbfolder) or die "ERROR: Cannot open $dbfolder: $!";
my @dbs = sort(grep(/^[0-9]+/, readdir(DBF)));

system("mv makebowtiedb.* 00.script/02.makebowtiedb.script/run.$run");
system("chmod 777 -R 00.script/02.makebowtiedb.script/run.$run");
system("chmod 777 -R $dbfolder");

system("rm -rf 00.script/$logfolder");
system("mkdir -p 00.script/$logfolder");

foreach my $db (@dbs){
	my $shell = "00.script/$logfolder/bowtie.$seqtype.$db.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N bowtie.$seqtype.$db.sh\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=40gb\n";
		print SHL "\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
		$thread = 2;
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	if($platform eq "sapelo"){
    	print SHL "module load bowtie2/2.2.4\n";
	}elsif($platform eq "zcluster"){
		print SHL "PATH=/usr/local/bowtie2/2.2.3/bin/:\$PATH\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	print SHL "mkdir -p $tgtfolder/$db\n";
			
	my @R1 = ();
	my @R2 = ();
	my @R3 = ();
	my @R4 = ();
	foreach my $sub (@subs){
		if($mode1 eq "paired-end" and $sub =~ /F$|Fu$|R$/){
			my $new = "$qryfolder/$sub/$sub.fna";
			push @R4, $new;
		}elsif($mode1 eq "paired-end"){
			my $new = "$qryfolder/$sub/$sub.R1.fasta_pairs_R1.fasta";
			push @R1, $new;
			$new = "$qryfolder/$sub/$sub.R2.fasta_pairs_R2.fasta";
			push @R2, $new;
			$new = "$qryfolder/$sub/$sub.R1.fasta_singles.fasta";
			push @R3, $new;
		}elsif($mode1 eq "single-end"){
			push @R1, "$qryfolder/$sub/$sub.fasta";
		}else{
			die "Error: Please specify the mode as 'single-end' or 'paired-end'.";
		}
	}
	my $unmap = "--no-unal";
	#if($run == 0){$unmap = "";}
	
	if($mode eq "end-to-end"){
		if($mode1 eq "paired-end"){
			print SHL "time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread -1 ", join(",", @R1), " -2 ", join(",", @R2), " -U ", join(",", @R3), " -S $tgtfolder/$db/bowtie.out.$db.sam\n";
		}elsif($mode1 eq "single-end"){
			print SHL "time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread -U ", join(",", @R1), " -S $tgtfolder/$db/bowtie.out.$db.sam\n";
		}else{
			die "Error: Please specify the mode as 'single-end' or 'paired-end'.";
		}
		if(scalar(@R4) > 0){
			print SHL "time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread -U ", join(",", @R4), " -S $tgtfolder/$db/bowtie.out.$db.long.sam\n";
		}
	}elsif($mode eq "local"){
		if($mode1 eq "paired-end"){
			print SHL "time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread --local -1 ", join(",", @R1), " -2 ", join(",", @R2), " -U ", join(",", @R3), " -S $tgtfolder/$db/bowtie.out.$db.sam\n";
		}elsif($mode1 eq "single-end"){
			print SHL "time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread --local -U ", join(",", @R1), " -S $tgtfolder/$db/bowtie.out.$db.sam\n";
		}else{
			die "Error: Please specify the mode as 'single-end' or 'paired-end'.";
		}
		if(scalar(@R4) > 0){
			print SHL "time bowtie2 -f -x $dbfolder/$db/$db $unmap -p $thread --local -U ", join(",", @R4), " -S $tgtfolder/$db/bowtie.out.$db.long.sam\n";
		}
	}
	
	print SHL "touch 00.script/$logfolder/$db.done.log\n";
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

close(QRY);
close(DBF);
system("echo 'Finished 03.bowtie.folder.pl!' >> job.monitor.txt");

system("chmod 777 -R 00.script/$logfolder");
