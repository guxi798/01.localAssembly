#!/usr/bin/perl -w
#running command: time perl 00.script/a9.folder_cut_read_by_length.pl ../03.Appalachian/01.data/02.Fasta 01.data/02.Fasta 50 fasta zcluster

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $read_len = shift @ARGV;
my $seq_format = lc(shift @ARGV);
my $platform = lc(shift @ARGV);
my $script = "a9.cut_length_script";

## start running the script
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^\w.+/, readdir(SRC)));
closedir SRC;

system("rm -rf 00.script/$script");
system("mkdir -p 00.script/$script");

foreach my $sub (@subs){
	my $shell = "00.script/$script/cut.read.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N cut.read.$sub.sh\n";
		print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=2gb\n";
		print SHL "\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
		
		print SHL "mkdir -p $tgtfolder/$sub\n";
		print SHL "time perl 00.script/a9.cut_read_by_length.pl $srcfolder/$sub/$sub.R1.fasta_pairs_R1.fasta $tgtfolder/$sub/$sub.R1.fasta_pairs_R1.fasta $read_len fasta\n";
		print SHL "time perl 00.script/a9.cut_read_by_length.pl $srcfolder/$sub/$sub.R2.fasta_pairs_R2.fasta $tgtfolder/$sub/$sub.R2.fasta_pairs_R2.fasta $read_len fasta\n";
		print SHL "time perl 00.script/a9.cut_read_by_length.pl $srcfolder/$sub/$sub.R1.fasta_singles.fasta $tgtfolder/$sub/$sub.R1.fasta_singles.fasta $read_len fasta\n";

		close(SHL);
		system("chmod u+x $shell");
		system("qsub $shell");
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
		
		print SHL "mkdir -p $tgtfolder/$sub\n";
		print SHL "time perl 00.script/a9.cut_read_by_length.pl $srcfolder/$sub/$sub.R1.fasta_pairs_R1.fasta $tgtfolder/$sub/$sub.R1.fasta_pairs_R1.fasta $read_len fasta\n";
		print SHL "time perl 00.script/a9.cut_read_by_length.pl $srcfolder/$sub/$sub.R2.fasta_pairs_R2.fasta $tgtfolder/$sub/$sub.R2.fasta_pairs_R2.fasta $read_len fasta\n";
		print SHL "time perl 00.script/a9.cut_read_by_length.pl $srcfolder/$sub/$sub.R1.fasta_singles.fasta $tgtfolder/$sub/$sub.R1.fasta_singles.fasta $read_len fasta\n";
		
		close(SHL);
		system("chmod u+x $shell");
		system("qsub -q rcc-30d $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}
