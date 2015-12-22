#!/usr/bin/perl -w
# run the script: time perl 00.script/07.blastx.back.pl 06.assembly/03.bowtie.nucl/run.0 01.data/05.splitGenes/01.Protein/run.0 07.map.back/03.bowtie.nucl/run.0 Sapelo 10

use strict;
system("echo 'Running 07.blastx.back.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $dbfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my ($run) = $srcfolder =~ /run\.([0-9]+)/;
my $thread = 1;

## start running the script
system("rm -rf 00.script/07.blastn.script/run.$run");
system("mkdir -p 00.script/07.blastn.script/run.$run");

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));

foreach my $sub (@subs){
	if($sub eq '99999'){next;}
	my $shell = "00.script/07.blastn.script/run.$run/blastn.back.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
    if($platform eq "sapelo"){
    	print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N blastn.back.$sub\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=1:00:00\n";
	    print SHL "#PBS -l mem=2gb\n";
    }elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	print SHL "\n";
    print SHL "cd \$PBS_O_WORKDIR\n";

	my $command1 = 0;
	if($platform eq "sapelo"){
		$command1 = "/usr/local/apps/ncbiblast+/2.2.29/bin";
    }elsif($platform eq "zcluster"){
		$command1 = "/usr/local/ncbiblast+/2.2.29/bin";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	print SHL "mkdir -p $tgtfolder/$sub\n";
	print SHL "time $command1/blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.out -evalue 1e-10  -outfmt 6 -num_threads 1 -max_target_seqs 1\n ";
	print SHL "time $command1/blastn -db $dbfolder/$sub/$sub.fasta -query $srcfolder/$sub/Trinity.new.fasta -out $tgtfolder/$sub/$sub.contigs.blast.xml.out -evalue 1e-10  -outfmt 5 -num_threads 1 -max_target_seqs 1\n ";
	
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

close SRC;
system("echo 'Finished 07.blastx.back.pl!' >> job.monitor.txt");
