#!/usr/bin/perl -w
# run the script: time perl ../../00.script/d14.blast.batch.pl tblastn protein.fa transcriptome.part3.v3.fa tblastn.batch 6 10000 zcluster

use strict;
use Bio::SeqIO;

my $mode = shift @ARGV;
my $srcfile = shift @ARGV;
my $dbfile = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $evalue = shift @ARGV;
my $size = shift @ARGV;
my $platform = lc(shift @ARGV);

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");

my $count = 0;
my $group = 1;
system("mkdir -p $tgtfolder/rh$group");
open(TGT, ">$tgtfolder/rh$group/rh$group.fasta");
print "Partition rh$group group.....\n";

while(my $seq_obj = $seqio_obj->next_seq){
	$count ++;
	
	if($count % $size == 0){
		my $shell = "rh$group.sh";
		open(SHL, ">$shell");
		if($platform eq "sapelo"){
			print SHL "#PBS -S /bin/bash\n";
			print SHL "#PBS -q batch\n";
			print SHL "#PBS -N blast.rh$group\n";
			print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
			print SHL "#PBS -l walltime=480:00:00\n";
			print SHL "#PBS -l mem=10gb\n\n";
			print SHL "cd \$PBS_O_WORKDIR\n\n";
			
			print SHL "time /usr/local/apps/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.out -outfmt 6 -num_threads 1 -max_target_seqs 10\n";
			print SHL "time /usr/local/apps/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.xml.out -outfmt 5 -num_threads 1 -max_target_seqs 10\n";

			system("chmod u+x $shell");
			system("qsub $shell");
		}elsif($platform eq "zcluster"){
			print SHL "#!/bin/bash\n\n";
			print SHL "time /usr/local/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.out -outfmt 6 -num_threads 1 -max_target_seqs 10\n";
			print SHL "time /usr/local/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.xml.out -outfmt 5 -num_threads 1 -max_target_seqs 10\n";

			system("chmod u+x $shell");
			system("qsub -q rcc-30d $shell");
		}else{
			die "Please provide the platform: Sapelo or Zcluster.\n";
		}

		$group ++;
		close TGT;
		system("mkdir -p $tgtfolder/rh$group");
		open(TGT, ">$tgtfolder/rh$group/rh$group.fasta");
		print "Partition rh$group group.....\n";
	}
	
	my $line = $seq_obj->id;
	my $desc = $seq_obj->desc;
	my $seq = $seq_obj->seq;
	
	print TGT ">$line $desc\n";
	print TGT "$seq\n";
	
}

my $shell = "rh$group.sh";
open(SHL, ">$shell");
if($platform eq "sapelo"){
	print SHL "#PBS -S /bin/bash\n";
	print SHL "#PBS -q batch\n";
	print SHL "#PBS -N blast.rh$group\n";
	print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
	print SHL "#PBS -l walltime=480:00:00\n";
	print SHL "#PBS -l mem=10gb\n\n";
	print SHL "cd \$PBS_O_WORKDIR\n\n";
	
	print SHL "time /usr/local/apps/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.out -outfmt 6 -num_threads 1 -max_target_seqs 10\n";
	print SHL "time /usr/local/apps/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.xml.out -outfmt 5 -num_threads 1 -max_target_seqs 10\n";

	system("chmod u+x $shell");
	system("qsub $shell");
}elsif($platform eq "zcluster"){
	print SHL "#!/bin/bash\n\n";
	print SHL "time /usr/local/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.out -outfmt 6 -num_threads 1 -max_target_seqs 10\n";
	print SHL "time /usr/local/ncbiblast+/2.2.29/bin/$mode -db $dbfile -query $tgtfolder/rh$group/rh$group.fasta -evalue $evalue -out $tgtfolder/rh$group/rh$group.xml.out -outfmt 5 -num_threads 1 -max_target_seqs 10\n";

	system("chmod u+x $shell");
	system("qsub -q rcc-30d $shell");
}else{
	die "Please provide the platform: Sapelo or Zcluster.\n";
}

$group ++;
close TGT;
system("mkdir -p $tgtfolder/rh$group");
open(TGT, ">$tgtfolder/rh$group/rh$group.fasta");
print "Partition rh$group group.....\n";


close TGT;
