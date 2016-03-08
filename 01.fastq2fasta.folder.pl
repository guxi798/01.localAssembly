#!/usr/bin/perl -w
# run the script: time perl 00.script/01.fastq2fasta.folder.pl 01.data/01.Fastq 01.data/02.Fasta Sapelo

use strict;
system("echo 'Running 01.fastq2fasta.folder.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $thread = 1;

## start running the script
opendir(SRC, $srcfolder) or die "Cannot open $srcfolder: $!";
my @subs = sort(grep(/\w+/, readdir(SRC)));

#system("mv makeblastdb.* 00.script/shell.script");
system("rm -rf 00.script/shell.script.previous");
system("mv 00.script/shell.script 00.script/shell.script.previous");
system("mkdir -p 00.script/shell.script");

foreach my $sub (@subs){
	my $shell = "00.script/shell.script/fastx.fastq2fasta.$sub.sh";
	open(SHL, ">$shell");
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N 00.script/shell.script/fastx.fastq2fasta.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=4:00:00\n";
	    print SHL "#PBS -l mem=20gb\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
    
    print SHL "\n";
    print SHL "cd \$PBS_O_WORKDIR\n";

	print SHL "mkdir -p $tgtfolder/$sub\n";
	if($sub =~ /F$|Fu$|R$/){
		print SHL "cp $srcfolder/$sub/$sub.fna $tgtfolder/$sub/\n";
	}else{
		print SHL "time awk '1 == (NR) % 4 || 2 == (NR) % 4' $srcfolder/$sub/$sub.R1.fastq_pairs_R1.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > $tgtfolder/$sub/$sub.R1.fasta_pairs_R1.fasta\n";
		print SHL "time awk '1 == (NR) % 4 || 2 == (NR) % 4' $srcfolder/$sub/$sub.R2.fastq_pairs_R2.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > $tgtfolder/$sub/$sub.R2.fasta_pairs_R2.fasta\n";
		print SHL "time awk '1 == (NR) % 4 || 2 == (NR) % 4' $srcfolder/$sub/$sub.R1.fastq_singles.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > $tgtfolder/$sub/$sub.R1.fasta_singles.fasta\n";
		print SHL "sed -i \"s/\\/1//g\" $tgtfolder/$sub/$sub.R1.fasta_singles.fasta\n";
		print SHL "sed -i \"s/\\/2//g\" $tgtfolder/$sub/$sub.R1.fasta_singles.fasta\n";
		#print SHL "time awk '1 == (NR) % 4 || 2 == (NR) % 4' $srcfolder/$sub/$sub.R1.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > $tgtfolder/$sub/$sub.R1.fasta\n";
		#print SHL "time awk '1 == (NR) % 4 || 2 == (NR) % 4' $srcfolder/$sub/$sub.R2.fastq | awk '{gsub(\"^@\", \">\", \$0); print \$0}' > $tgtfolder/$sub/$sub.R2.fasta\n";
	}
	
	close(SHL);
		
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "scluster"){
		system("qsub -q rcc-30d -pe thread $thread $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

close(SRC);
system("echo 'Finished 01.fastq2fasta.folder.pl!' >> job.monitor.txt");
system("chmod 777 -R 01.data/01.Fastq");
system("chmod 777 -R 01.data/02.Fasta");
system("chmod 777 -R 00.script");
