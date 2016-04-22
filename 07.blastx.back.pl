#!/usr/bin/perl -w
# run the script: time perl 00.script/07.blastx.back.pl 06.assembly/03.bowtie.nucl/run.0 01.data/05.splitGenes/01.Protein/run.0 07.map.back/03.bowtie.nucl/run.0 Sapelo 10

use strict;
system("echo 'Running 07.blastx.back.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $dbfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $srcfolder =~ /run\.([0-9]+)/;
my $thread = 1;

#=pod
## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
#if($run == run.0){
#	push @chks, 99999;
#}

while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("truncate.header.$chk.o*");
		my @stdout = glob("truncate.header.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' truncate.header.$chk.e* >> 00.script/06.truncate.script/run.$run/summary.error.log\n");
			#system("echo 'success' > 00.script/shell.script/truncate.header.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/06.truncate.script/run.$run/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/Trinity.new.fasta")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					system("rm -f truncate.header.$chk.*");
					system("echo 'Resubmitting the job: truncate.header.$chk.sh' >> job.monitor.txt");
					system("qsub 00.script/06.truncate.script/run.$run/truncate.header.$chk.sh");
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
#=cut

## start running the script
system("mv truncate.header.* 00.script/06.truncate.script/run.$run/");
system("chmod 777 -R 00.script/06.truncate.script/run.$run");
system("chmod 777 -R $srcfolder");

system("rm -rf 00.script/07.blastx.script/run.$run");
system("mkdir -p 00.script/07.blastx.script/run.$run");

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));

foreach my $sub (@subs){
	if($sub eq '99999'){next;}
	my $shell = "00.script/07.blastx.script/run.$run/blastx.back.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
    if($platform eq "sapelo"){
    	print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N blastx.back.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=1:00:00\n";
	    print SHL "#PBS -l mem=2gb\n";
		print SHL "\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
   }elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	my $command1 = 0;
	if($platform eq "sapelo"){
		$command1 = "/usr/local/apps/ncbiblast+/2.2.29/bin";
    }elsif($platform eq "zcluster"){
		$command1 = "/usr/local/ncbiblast+/2.2.29/bin";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	print SHL "mkdir -p $tgtfolder/$sub\n";
	print SHL "time $command1/blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-5  -outfmt 6 -num_threads 1 -max_target_seqs 1\n ";
	print SHL "time $command1/blastx -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-5  -outfmt 5 -num_threads 1 -max_target_seqs 1\n ";

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

close SRC;
system("echo 'Finished 07.blastx.back.pl!' >> job.monitor.txt");

system("chmod 777 -R 00.script/07.blastx.script/run.$run");