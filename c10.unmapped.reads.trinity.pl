#!/usr/bin/perl -w
# run the script: time perl 00.script/c10.unmapped.reads.trinity.pl 09.bowtie.full.length $platform

use strict;

my $srcfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $thread = 16;
my $memory = 160;

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^\w+/, readdir(SRC)));
closedir SRC;
my @folder = map "", 1..scalar(@subs);
my @left = map {"$srcfolder/".$_."/unmapped.reads.".$_.".1.fasta"} @subs;
my @right = map {"$srcfolder/".$_."/unmapped.reads.".$_.".2.fasta"} @subs;

my $shell = "00.script/c10.unmapped.reads.trinity.sh";
        open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
if($platform eq "sapelo"){
    print SHL "#PBS -S /bin/bash\n";
    print SHL "#PBS -q batch\n";
    print SHL "#PBS -N c10.unmapped.reads.trinity\n";
    print SHL "#PBS -l nodes=1:ppn=$thread:HIGHMEM\n";
    print SHL "#PBS -l walltime=120:00:00\n";
    print SHL "#PBS -l mem=",$memory,"gb\n\n";

    print SHL "cd \$PBS_O_WORKDIR\n";
    print SHL "module load trinity/r20140717\n";

    print SHL "time /usr/local/apps/trinity/r20140717/Trinity --seqType fa --CPU $thread --JM ",$memory,"G --left ", join(",",@left), " --right ", join(",",@right), " --output 10.unmapped.reads.trinity\n";
	close(SHL);

    system("chmod u+x $shell");
    system("qsub $shell");
}elsif($platform eq "zcluster"){
	print SHL "#!/bin/bash\n";
	print SHL "export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:\${LD_LIBRARY_PATH}\n";
	print SHL "export PATH=/usr/local/gmap-gsnap/latest/bin/:\${PATH}\n\n";
	$thread = 4;
	$memory = "80";
	
	print SHL "time /usr/local/trinity/r20140717/Trinity --seqType fa --CPU $thread --JM ",$memory,"G --left ", join(",",@left), " --right ", join(",",@right), " --output 10.unmapped.reads.trinity\n";

	close SHL;
	system("chmod u+x $shell");
	system("qsub -q rcc-m128-30d -pe thread $thread -l mem_total=$memory"."g $shell");
}else{
	die "Please provide the platform: 'Sapelo' or 'Zcluster'";
}
