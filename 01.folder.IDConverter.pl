#!/usr/bin/perl -w
# run the script: time perl 00.script/01.folder.IDConverter.pl  01.data/02.Fasta Sapelo

use strict;
system("echo 'Running 01.folder.IDConverter.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $thread = 1;

## start running the script
system("mv fastaCombinePairedEnd.* 00.script/shell.script");
system("rm -rf 00.script/shell.script.previous");
system("mv 00.script/shell.script 00.script/shell.script.previous");
system("mkdir -p 00.script/shell.script");

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/\w+/, readdir(SRC)));

foreach my $sub (@subs){
	my $shell = "00.script/shell.script/IDCoverter.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
    print SHL "#PBS -S /bin/bash\n";
    print SHL "#PBS -q batch\n";
    print SHL "#PBS -N IDCoverter.$sub\n";
    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
    print SHL "#PBS -l walltime=1:00:00\n";
    print SHL "#PBS -l pvmem=20gb\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
   
    print SHL "\n";
    print SHL "cd \$PBS_O_WORKDIR\n";
	if($sub =~ /F$|Fu$|R$/){
    	print SHL "time cat $srcfolder/$sub/$sub.fna | awk '{if(\$_ ~ /^>/){count++; print \">\"count} else{print \$_}}' > $srcfolder/$sub/$sub.simple.fasta\n";
    }else{
		print SHL "time cat $srcfolder/$sub/$sub.R1.fasta_pairs_R1.fasta | awk '{if(NR%2==0){print \$_} else{print \">\"NR/2+0.5}}' > $srcfolder/$sub/$sub.R1.fasta_simple.fasta\n";
		print SHL "time cat $srcfolder/$sub/$sub.R2.fasta_pairs_R2.fasta | awk '{if(NR%2==0){print \$_} else{print \">\"NR/2+0.5}}' > $srcfolder/$sub/$sub.R2.fasta_simple.fasta\n";
		print SHL "time cat $srcfolder/$sub/$sub.R1.fasta_singles.fasta | awk '{if(NR%2==0){print \$_} else{print \">\"NR/2+0.5}}' > $srcfolder/$sub/$sub.singles.fasta_simple.fasta\n";
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

system("echo 'Finished 01.folder.IDConverter.pl!' >> job.monitor.txt");
