#!/usr/bin/perl -w
# run the script: time perl 00.script/10.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.0 07.map.back/03.bowtie.nucl/run.0 01.data/05.splitGenes/02.Transcript/run.1 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 

use strict;
use Bio::SeqIO;

system("echo 'Running 10.transfer.saturate.seq.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;			## assembled contig folder
my $blastfolder = shift @ARGV;			## blast result folder
my $tgtfolder = shift @ARGV;			## target folder to put contigs for new run
my $reffile = shift @ARGV;				## reference file with gene ID and protein length info
my $mode = shift @ARGV;					## abs: absolute value; pct: percent value
my $cutoff = shift @ARGV;				## absolute AA number, or percent
my $sleeptime = shift @ARGV;
my ($run) = $srcfolder =~ /\/run\.([0-9]+)/;
my $errfile = "00.script/10.transfer.script/run.$run/transfer.saturate.seq.e";
my $outfile = "00.script/10.transfer.script/run.$run/transfer.saturate.seq.o";

#=pod
## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("blastx.back.$chk.o*");
		my @stdout = glob("blastx.back.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' blastx.back.$chk.e* >> 00.script/07.blastx.script/run.$run/summary.error.log\n");
			#system("echo 'success' > 00.script/shell.script/blastx.back.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/07.blastx.script/run.$run/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$blastfolder/$chk/$chk.contigs.blast.out")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					system("rm -f blastx.back.$chk.*");
					system("echo 'Resubmitting the job: truncate.header.$chk.sh' >> job.monitor.txt");
					system("qsub 00.script/07.blastx.script/run.$run/blastx.back.$chk.sh");
					#last; # There are some jobs failed
				}
				else{
					$i++;
				}
			}
			if($i == scalar @chks){
				system("echo 'All jobs have been finished successfully' >> job.monitor.txt"); # error file is empty
				last;
			}
		}else{
			die "ERROR: something went wrong in previous steps\n";	
		}
	}
	sleep $sleeptime;
}
close CHK;

#=cut

## start running the script
system("mv blastx.back.* 00.script/07.blastx.script/run.$run/");
system("mv blastn.back.* 00.script/07.blastn.script/run.$run/");
system("chmod 777 -R 00.script/07.blastx.script/run.$run");
system("chmod 777 -R 00.script/07.blastn.script/run.$run");
system("chmod 777 -R $blastfolder");
system("chmod 777 -R 07.map.back/02.blastn/run.$run");

system("rm -rf 00.script/10.transfer.script/run.$run");
system("mkdir -p 00.script/10.transfer.script/run.$run");
system("mkdir -p $tgtfolder");

opendir(SRC, $blastfolder) or die "ERROR: Cannot open $blastfolder: $!";
my @subs = sort(grep(/^[0-9]+$/, readdir SRC));

open(ERR, ">$errfile") or die "ERROR: Cannot write $errfile: $!";
open(OUT, ">$outfile") or die "ERROR: Cannot write $outfile: $!";
foreach my $sub (@subs){
	print OUT "$sub\n";
	open(BLS, "$blastfolder/$sub/$sub.contigs.blast.out") or die "ERROR: Cannot open $blastfolder/$sub/$sub.contigs.blast.out: $!";
	my %valid = ();
	foreach my $line (<BLS>){
		chomp $line;
		my @lines = split(/\t/, $line);
		my ($contig, $gene) = @lines[0..1];
		
		if(not exists $valid{$contig}){
            $valid{$contig} = 0;
        }
		
	}
	close(BLS);
	
	open(TRI, "$srcfolder/$sub/Trinity.new.fasta") or die "ERROR: Cannot open $srcfolder/$sub/Trinity.new.fasta: $!";
	system("mkdir -p $tgtfolder/$sub");
	open(TGT2, ">$tgtfolder/$sub/$sub.fasta") or die "ERROR: Cannot write $tgtfolder/$sub/$sub.fasta: $!";
	my $mark = 0;
	foreach my $line (<TRI>){
		chomp $line;
		my @lines = split(/\s+/, $line);
		if($lines[0] =~ /^>/){
			$lines[0] =~ s/>//;
		if(exists $valid{$lines[0]}){
				print TGT2 ">$lines[0] $lines[1]\n";
				$mark =2;
			}else{
				$mark =0;
			}
		}else{
			if($mark == 2){
				print TGT2 "$line\n";
			}
		}
	}
}

close TRI;
close TGT2;
close ERR;
close OUT;
close SRC;

system("grep -E 'ERROR|Error|error' 00.script/10.transfer.script/run.$run/transfer.saturate.seq.e > 00.script/10.transfer.script/run.$run/summary.error.log");
system("echo 'success' > 00.script/10.transfer.script/run.$run/transfer.saturate.seq.log");

system("echo 'Finished 10.transfer.saturate.seq.pl!' >> job.monitor.txt");

system("chmod 777 -R 00.script/10.transfer.script/run.$run");
