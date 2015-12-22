#!/usr/bin/perl -w
# run the script: time perl 00.script/040.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.0 04.retrieve.reads/03.bowtie.nucl/run.0 nucl bowtie.log/bowtie.run.0 1000 genome Sapelo 10

use strict;
system("echo 'Running 040.folder.retrievebowtie.reads.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $seqtype = shift @ARGV;
my $logfolder = shift @ARGV;
my $blocksize = shift @ARGV;
my $scale = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $srcfolder =~ /run\.([0-9]+)/;
my $thread = 1;

## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
@chks = @chks;

opendir(QRY, "01.data/02.Fasta") or die "ERROR: Cannot open 01.data/02.Fasta: $!";
my @sams = sort(grep(/^\w.+/, readdir(QRY)));
@sams = @sams;

while(1){
	my $count = 0;
	my @temp = @chks;
	my @temp2 = @sams;
	my $i = 0;
	while(my $chk = shift @temp){
		@temp2 = @sams;
		while(my $sam = shift @temp2){
			my @stderr = glob("bowtie.$seqtype.$chk.$sam.o*");
			my @stdout = glob("bowtie.$seqtype.$chk.$sam.e*");
			if(!($stderr[0] and $stdout[0])){
				last; # job hasn't finished, no log file
			}else{
				system("grep -E 'ERROR|Error|error' bowtie.$seqtype.$chk.$sam.e* >> 00.script/$logfolder/summary.error.log\n");
				$count ++; # job has finished, add one count
			}
		}			
	}
	@temp = @chks;
	if($count == (scalar @chks) * (scalar @sams)){ # all jobs have finished
		if(!-s "00.script/$logfolder/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
				system("echo 'All jobs have been finished successfully' >> job.monitor.txt"); # error file is empty
				last;
		}else{
			die "ERROR: something went wrong in previous steps\n";	
		}
	}
	sleep $sleeptime;
}
close CHK;
close QRY;

## start running the script

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));
#if($run eq '0' and $scale eq 'genome'){	# assemble unmapped reads in the first run
#	push @subs, 99999;
#}

system("mv bowtie.* 00.script/$logfolder");
system("rm -rf 00.script/04.retrieve.script/run.$run");
system("mkdir -p 00.script/04.retrieve.script/run.$run");
system("rm -rf $tgtfolder");

foreach my $sub (@subs){
	my $shell = "00.script/04.retrieve.script/run.$run/retrievebowtie.reads.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N retrievebowtie.reads.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=12:00:00\n";
	    print SHL "#PBS -l mem=30gb\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

    print SHL "\n";
    print SHL "cd \$PBS_O_WORKDIR\n";
	print SHL "mkdir -p $tgtfolder/$sub\n";
	
	my @tempsams = @sams;
	if($sub ne "99999"){
		foreach my $sam (@tempsams){
			if($sam =~ /F$|Fu$|R$/) { ## 454 data
				print SHL "time perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.tab 01.data/02.Fasta/$sam/$sam.simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.long.fasta\n";
			}else{
				print SHL "time perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R1.tab 01.data/02.Fasta/$sam/$sam.R1.fasta_simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta\n";
				print SHL "time perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.R2.tab 01.data/02.Fasta/$sam/$sam.R2.fasta_simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R2.fasta\n";
				print SHL "time perl 00.script/040.retrievebowtie.reads.pl $srcfolder/$sub/bowtie.out.$sub.$sam.single.tab 01.data/02.Fasta/$sam/$sam.singles.fasta_simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta\n";
			}
		}
	}else{
		foreach my $sam (@tempsams){
			if($sam =~ /F$|Fu$|R$/){
				print SHL "time perl 00.script/040.retrieveunmap.reads.pl $srcfolder $sam F 01.data/02.Fasta/$sam/$sam.simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.long.fasta\n";
			}else{
				print SHL "time perl 00.script/040.retrieveunmap.reads.pl $srcfolder $sam R1 01.data/02.Fasta/$sam/$sam.R1.fasta_simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta\n";
				print SHL "time perl 00.script/040.retrieveunmap.reads.pl $srcfolder $sam R2 01.data/02.Fasta/$sam/$sam.R2.fasta_simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R2.fasta\n";
				print SHL "time perl 00.script/040.retrieveunmap.reads.pl $srcfolder $sam S 01.data/02.Fasta/$sam/$sam.singles.fasta_simple.fasta $blocksize >> $tgtfolder/$sub/retrieved.$sub.R1.fasta\n";
			}
		}
	}
	print SHL "sed -i 's/\\(>[A-Za-z0-9_]*\\)/\\1\\/1/' $tgtfolder/$sub/retrieved.$sub.R1.fasta\n";
	print SHL "sed -i 's/\\(>[A-Za-z0-9_]*\\)/\\1\\/2/' $tgtfolder/$sub/retrieved.$sub.R2.fasta\n";
	
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
system("echo 'Finished 040.folder.retrievebowtie.reads.pl!' >> job.monitor.txt");
