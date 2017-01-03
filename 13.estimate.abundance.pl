#!/usr/bin/perl -w
# run the script: time perl 00.script/13.estimate.abundance.pl 01.data/01.Fastq 13.abundance/run.3 01.data/05.splitGenes/02.Transcript/run.3/contigs.good.fasta single-end

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;			## reads file
my $tgtfolder = shift @ARGV;			## target folder to put contigs for new run
my $reffile = shift @ARGV;				## reference file with gene ID and protein length info
my $method = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $trinity = shift @ARGV;
my $thread = 4;
if(!$trinity){$trinity = "";}

## read in blast result file
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^\w+/, readdir SRC));
closedir SRC;
print "$srcfolder\n";
print join("\n", @subs), "\n";

system("rm -rf 00.script/13.rsem.script");
system("mkdir -p 00.script/13.rsem.script");

foreach my $sub (@subs){	## loop over parallized groups
    if($sub =~ /F$|Fu$|R$/){next;}
	my $shell = "00.script/13.rsem.script/rsem.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N 00.script/shell.script/rsem.bowtie2.single.$sub\n";
		print SHL "#PBS -l nodes=1:ppn=2:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=15gb\n";
		print SHL "\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

    my $command1 = 0;
	my $t = $thread / 2;
	if($platform eq "sapelo"){
		$command1 = "/usr/local/apps/trinity/r20140717";
	    print SHL "module load trinity/r20140717\n\n";
    }elsif($platform eq "zcluster"){
		$command1 = "/usr/local/trinity/r20140717";
		print SHL "export LD_LIBRARY_PATH=/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:\${LD_LIBRARY_PATH}\n";
		print SHL "export PATH=/usr/local/gmap-gsnap/latest/bin/:\${PATH}\n";
		print SHL "export PATH=\${PATH}:/usr/local/cd-hit/4.6.1-2012-08-27/:/usr/local/rsem/1.2.20/:/usr/local/express/1.5.1/\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	if($mode eq "paired-end"){
		print SHL "time $command1/util/align_and_estimate_abundance.pl --transcripts $reffile --seqType fq --left $srcfolder/$sub/$sub.R1.fastq_pairs_R1.fastq --right $srcfolder/$sub/$sub.R2.fastq_pairs_R2.fastq --single $srcfolder/$sub/$sub.R1.fastq_singles.fastq --output_dir $tgtfolder/$sub --est_method $method --aln_method bowtie2 $trinity\n";
	}elsif($mode eq "single-end"){
		print SHL "time $command1/util/align_and_estimate_abundance.pl --transcripts $reffile --seqType fq --single $srcfolder/$sub/$sub.fastq --output_dir $tgtfolder/$sub --est_method $method --aln_method bowtie2 $trinity\n";
	}
	
	close SHL;
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "zcluster"){
		system("qsub -q rcc-30d -pe thread $t $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

}

