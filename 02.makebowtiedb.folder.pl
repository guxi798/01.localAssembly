#!/usr/bin/perl -w
# run the script: time perl 00.script/02.makebowtiedb.folder.pl 01.data/05.splitGenes/02.Transcript/run.0 Sapelo 10

use strict;
system("echo 'Running 02.makebowtiedb.folder.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $thread = 1;
my ($run) = $srcfolder =~ /run\.([0-9]+)/;

## start running the script
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));

system("rm -rf 00.script/02.makebowtiedb.script/run.$run");
system("mkdir -p 00.script/02.makebowtiedb.script/run.$run");

foreach my $sub (@subs){
	my $shell = "00.script/02.makebowtiedb.script/run.$run/makebowtiedb.$sub.sh";
	my $base = $sub;
	$base =~ s/\.fasta//;
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N makebowtiedb.$sub.sh\n";
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
		print SHL "PATH=/usr/local/bowtie2/2.2.3/bin/:\$PATH\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	print SHL "cd $srcfolder/$sub\n";
	print SHL "time bowtie2-build -f -q $sub.fasta $base\n";
	print SHL "cd ../../../../../\n";
	print SHL "touch 00.script/02.makebowtiedb.script/run.$run/$sub.done.log\n";
		
	close(SHL);
		
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "zcluster"){
		system("qsub -q rcc-30d $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

close(SRC);
system("echo 'Finished 02.makebowtiedb.folder.pl!' >> job.monitor.txt");
system("chmod 777 -R 00.script/02.makebowtiedb.script/");
system("chmod 777 -R 01.data/05.splitGenes");