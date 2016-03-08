#!/usr/bin/perl -w
# run the script: time perl 00.script/13.estimate.abundance.pl 01.data/01.Fastq 13.abundance/run.3 01.data/05.splitGenes/02.Transcript/run.3/contigs.good.fasta

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;			## reads file
my $tgtfolder = shift @ARGV;			## target folder to put contigs for new run
my $reffile = shift @ARGV;				## reference file with gene ID and protein length info

## read in blast result file
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^\w+$/, readdir SRC));

system("rm -rf 00.script/shell.script.rsem");
system("mkdir -p 00.script/shell.script.rsem");

foreach my $sub (@subs){	## loop over parallized groups
    if($sub =~ /F$|Fu$|R$/){next;}
	my $shell = "00.script/shell.script.rsem/rsem.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
    print SHL "#PBS -S /bin/bash\n";
    print SHL "#PBS -q batch\n";
    print SHL "#PBS -N 00.script/shell.script/rsem.bowtie2.single.$sub\n";
    print SHL "#PBS -l nodes=1:ppn=2:AMD\n";
    print SHL "#PBS -l walltime=48:00:00\n";
    print SHL "#PBS -l mem=15gb\n";

    print SHL "\n";
    print SHL "module load trinity/r20140717\n\n";
    print SHL "cd \$PBS_O_WORKDIR\n";

    print SHL "#!/bin/bash\n";
	print SHL "time /usr/local/apps/trinity/r20140717/util/align_and_estimate_abundance.pl --transcripts $reffile --seqType fq --left $srcfolder/$sub/$sub.R1.fastq_pairs_R1.fastq --right $srcfolder/$sub/$sub.R2.fastq_pairs_R2.fastq --output_dir $tgtfolder/$sub --est_method RSEM --aln_method bowtie2 --trinity_mode\n";
		
	#print SHL "time /usr/local/apps/trinity/r20140717/util/align_and_estimate_abundance.pl --transcripts $reffile --seqType fq --single $srcfolder/$sub/$sub.single --output_dir $tgtfolder/$sub --est_method RSEM --aln_method bowtie2 --trinity_mode\n";

	close SHL;
	system("chmod u+x $shell");
	system("qsub $shell");

}

