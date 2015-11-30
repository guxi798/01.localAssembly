#!/usr/bin/perl -w
# run the script: time perl 00.script/c1.trinity.pl

use strict;

my $srcfolder = "01.data/01.Fastq";

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^\w+/, readdir(SRC)));
my @folder = map "", 1..scalar(@subs);
my @left = map {"01.data/01.Fastq/".$_."/".$_.".R1.ncRNA"} @subs;
my @right = map {"01.data/01.Fastq/".$_."/".$_.".R2.ncRNA"} @subs;

my $shell = "00.script/c1.trinity.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
    print SHL "#PBS -S /bin/bash\n";
    print SHL "#PBS -q batch\n";
    print SHL "#PBS -N 00.script/trinity_sub\n";
    print SHL "#PBS -l nodes=n26:ppn=30:HIGHMEM\n";
    print SHL "#PBS -l walltime=480:00:00\n";
    print SHL "#PBS -l mem=800gb\n";
    	
    print SHL "#PBS -M guxi798\@hotmail.com\n";
   	print SHL "#PBS -m abe\n\n";
     
    print SHL "module load trinity/r20140717\n";
    print SHL "cd \$PBS_O_WORKDIR\n";
	print SHL "time /usr/local/apps/trinity/r20140717/Trinity --seqType fq --left ", join(",",@left), " --right ", join(",",@right), " --normalize_reads --CPU 30 --JM 800G --output 02.trinity\n";
	
	close(SHL);
	
	system("chmod u+x $shell");
	system("qsub $shell");

close(SRC);
