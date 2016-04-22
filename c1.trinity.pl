#!/usr/bin/perl -w
# run the script: time perl 00.script/c1.trinity.pl

use strict;

my $srcfolder = "01.data/02.Fasta";
my $platform = lc(shift @ARGV);
my $thread = 8;
my $memory = 80;

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^\w+/, readdir(SRC)));
my @folder = map "", 1..scalar(@subs);
my @left = map {"$srcfolder/".$_."/".$_.".R1.fasta"} @subs;
my @right = map {"$srcfolder/".$_."/".$_.".R2.fasta"} @subs;

my $shell = "00.script/c1.trinity.sh";
open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";

if($platform eq "sapelo"){
    print SHL "#PBS -S /bin/bash\n";
    print SHL "#PBS -q batch\n";
    print SHL "#PBS -N 00.script/trinity_sub\n";
    print SHL "#PBS -l nodes=n26:ppn=$thread:HIGHMEM\n";
    print SHL "#PBS -l walltime=480:00:00\n";
    print SHL "#PBS -l mem=",$memory,"gb\n";
    	
    print SHL "#PBS -M guxi798\@hotmail.com\n";
   	print SHL "#PBS -m abe\n\n";
     
    print SHL "module load trinity/r20140717\n";
    print SHL "cd \$PBS_O_WORKDIR\n";
	print SHL "time /usr/local/apps/trinity/r20140717/Trinity --seqType fa --left ", join(",",@left), " --right ", join(",",@right), " --normalize_reads --CPU $thread --JM ", $memory, "G --output 05.Trinity\n";
	close(SHL);
	
	system("chmod u+x $shell");
	system("qsub $shell");
}elsif($platform eq "zcluster"){
	print SHL "#!/bin/bash\n\n";
	print SHL "export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:\${LD_LIBRARY_PATH}\n";
	print SHL "export PATH=/usr/local/gmap-gsnap/latest/bin/:\${PATH}\n\n";
	$thread = 4;
	$memory = "50";
	print SHL "time /usr/local/trinity/r20140717/Trinity --seqType fa --left ", join(",",@left), " --right ", join(",",@right), " --normalize_reads --CPU $thread --JM ", $memory, "G --output 05.Trinity\n";
	
	close(SHL);
	
	system("chmod u+x $shell");
	system("qsub -q rcc-30d -pe thread $thread -l mem_total=$memory"."g $shell");
}
	

close(SRC);
