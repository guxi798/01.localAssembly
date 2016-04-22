#!/usr/bin/perl -w
#running command: time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 08.full.length/full.length.contigs 09.unmapped.reads Sapelo

use strict;

## read in parameters required by the script
my $qryfolder = shift @ARGV;
my $dbfile = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $thread = 4;

## start running the script
opendir(QRY, $qryfolder) or die "ERROR: Cannot open $qryfolder: $!";
my @subs = sort(grep(/^\w.+/, readdir(QRY)));
closedir QRY;

system("rm -rf 00.script/bowtie.log/bowtie.run.c10");
system("mkdir -p 00.script/bowtie.log/bowtie.run.c10");

foreach my $sub (@subs){
	my $shell = "00.script/bowtie.log/bowtie.run.c10/bowtie.full.length.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N bowtie.full.length.$sub.sh\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=40gb\n";
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
		print SHL "PATH=/usr/local/bowtie2/2.2.3/bin/:\$PATH";
		$thread /= 2;
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	print SHL "mkdir -p $tgtfolder/$sub\n";
			
	my $R1 = "$qryfolder/$sub/$sub.R1.fasta_pairs_R1.fasta";
	my $R2 = "$qryfolder/$sub/$sub.R2.fasta_pairs_R2.fasta";
	my $R3 = "$qryfolder/$sub/$sub.R1.fasta_singles.fasta";
	my $R4 = "$qryfolder/$sub/$sub.fna";
	
	if($sub =~ /F$|Fu$|R$/){
		print SHL "time bowtie2 -f -x $dbfile -p $thread -U $R4 --un $tgtfolder/$sub/unmapped.reads.$sub.long.fasta\n";
	}else{
		print SHL "time bowtie2 -f -x $dbfile -p $thread -1 $R1 -2 $R2 -U $R3 --un-conc $tgtfolder/$sub/unmapped.reads.$sub.fasta --un $tgtfolder/$sub/unmapped.reads.$sub.single.fasta\n";
		print SHL "cat $tgtfolder/$sub/unmapped.reads.$sub.single.fasta >> $tgtfolder/$sub/unmapped.reads.$sub.1.fasta";
	}
	
	print SHL "\n";
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
