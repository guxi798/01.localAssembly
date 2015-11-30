#!/usr/bin/perl -w
# run the script: time perl 00.script/13.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.2/01.trinity 07.map.back/03.bowtie.nucl/run.2 01.data/05.splitGenes/02.Transcript/run.100 01.data/04.GeneOfInterest/GeneID.v1.txt pct 1 10 

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
my $errfile = "00.script/shell.script/transfer.saturate.seq.e";
my $outfile = "00.script/shell.script/transfer.saturate.seq.o";
my ($run) = $srcfolder =~ /\/(run\.[0-9]+)\//;

## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/02.Transcript/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("plot.blastx.$chk.o*");
		my @stdout = glob("plot.blastx.$chk.e*");
		if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			#system("grep -E 'ERROR|Error|error' plot.blastx.$chk.e* >> 00.script/shell.script/summary.error.log\n");
			system("touch 00.script/shell.script/summary.error.log");
			#system("echo 'success' > 00.script/shell.script/plot.blastx.$chk.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "00.script/shell.script/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "09.plot/03.bowtie.nucl/$run/$chk/$chk.evalue.emf" and -s "09.plot/03.bowtie.nucl/$run/$chk/$chk.identity.emf" and -s "09.plot/03.bowtie.nucl/$run/$chk/$chk.length.emf")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f plot.blastx.$chk.*");
					system("echo 'Resubmitting the job: plot.blastx.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/shell.script/plot.blastx.$chk.sh");
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

## start running the script
## prepare file and folder
system("mv plot.blastx.* 00.script/shell.script/");
system("mkdir -p $tgtfolder");
system("rm -rf 00.script/shell.script.previous");
system("mv 00.script/shell.script 00.script/shell.script.previous");
system("mkdir -p 00.script/shell.script");
## record file: error and output file
open(ERR, ">$errfile") or die "ERROR: Cannot write $errfile: $!";
open(OUT, ">$outfile") or die "ERROR: Cannot write $outfile: $!";

open(TGT1, ">$tgtfolder/full.length.contigs.fasta") or die "ERROR: Cannot write $tgtfolder/full.length.contigs.fasta: $!";

## read ref file and put gene ID and protein length into hash
open(REF, $reffile) or die "ERROR: Cannot open $reffile: $!";
<REF>; #remove head line

my %hash = ();
foreach my $line (<REF>){
	chomp $line;
	my @lines = split(/\t/, $line);
	if(not exists $hash{$lines[0]}){
		$hash{$lines[0]} = $lines[2];		# key: gene ID, value: proten len
	}
}

## read in blast result file
opendir(SRC, $blastfolder) or die "ERROR: Cannot open $blastfolder: $!";
my @subs = sort(grep(/^[0-9]+$/, readdir SRC));

foreach my $sub (@subs){	## loop over parallized groups
	print OUT "$sub\n";
	open(BLS, "$blastfolder/$sub/$sub.contigs.blast.out") or die "ERROR: Cannot open $blastfolder/$sub/$sub.contigs.blast.out: $!";
	my %full = ();
	foreach my $line (<BLS>){		## loop over blast record
		chomp $line;
		my @lines = split(/\t/, $line);
		my ($contig_info, $gene) = @lines[0..1];
		my @contig_info = split(/\|/, $contig_info);
		my $contig = $contig_info[0];
		my $left = ($lines[8], $lines[9])[$lines[8] > $lines[9]];		## subject start
		my $right = ($lines[8], $lines[9])[$lines[8] < $lines[9]];		## subject end
		
		if(not exists $full{$gene}){
			$full{$gene} = [0, abs($right-$left)+1, $hash{$gene}, $contig];
			#print abs($right-$left)+1, ", $hash{$gene}\n";
			if($left == 1 and $right == $hash{$gene}){
				${$full{$gene}}[0] = 1;
				#print "$contig, $left, $right, $hash{$gene}\n";
			}
		}else{
			push @{$full{$gene}}, $contig;			## there are multiple contigs mapped to the same gene
			if($left == 1 and $right == $hash{$gene}){
				${$full{$gene}}[0] = 1;
				${$full{$gene}}[1] = abs($right-$left)+1;
				${$full{$gene}}[2] = $hash{$gene};
			}
		}
	}
	close(BLS);
	
	my %contigs = ();
	my %incomplete = ();
	
	foreach my $key (sort(keys %full)){
		my $mark = shift @{$full{$key}}; 
		my $qlen = shift @{$full{$key}};	## query length
		my $slen = shift @{$full{$key}};	## subject length

		## the contig covers full length protein
		if($mark == 1){			
			#print "$mark, $qlen, $slen, ", join(", ", @{$full{$key}}),"\n";
			## remove all contigs under the same gene, not give chance to alt contig to grow 
			foreach my $contig (@{$full{$key}}){
				if(not exists $contigs{$contig}){
					$contigs{$contig} = [$key, $qlen, $slen, $mark];
				}elsif($key eq ${$contigs{$contig}}[0]){
					print ERR "Warnings: $contig is mapped to the same gene with multiple hit.\n";
				}else{
					print ERR "ERROR: something unexpected happened for $contig\n";
				}
			}
		}else{	## the contig covers only partial length protein
			foreach my $contig (@{$full{$key}}){
				if(not exists $incomplete{$contig}){
					$incomplete{$contig} = [$key, $qlen, $slen, $mark];
				}elsif($key eq ${$incomplete{$contig}}[0]){
					print ERR "Warnings: $contig is mapped to the same gene with multiple hit.\n";
				}else{
					print ERR "ERROR: something unexpected happened for $contig\n";
				}
			}
		}
		## the contig covers full length protein
	}
	
	#open(TRI, "$srcfolder/$sub/Trinity.new.fasta") or die "ERROR: Cannot open $srcfolder/$sub/Trinity.new.fasta: $!";
	system("mkdir -p $tgtfolder/$sub");
	open(TGT2, ">$tgtfolder/$sub/$sub.fasta") or die "ERROR: Cannot write $tgtfolder/$sub/$sub.fasta: $!";
	my $seqfile = "$srcfolder/$sub/Trinity.new.fasta";
	#my $outfile1 = "$tgtfolder/full.length.contigs.fasta";
	#my $outfile2 = "$tgtfolder/$sub/$sub.fasta";
	
	my $seqio_obj = Bio::SeqIO->new(-file => $seqfile, -format => "fasta");
	#my $seqout1 = Bio::SeqIO->new(-file => ">$outfile1", -format => "fasta");
	#my $seqout2 = Bio::SeqIO->new(-file => ">$outfile2", -format => "fasta");
	
	while(my $seq_obj = $seqio_obj->next_seq){
		my $line = $seq_obj->id;
		#print "$line\n";
		my @lines = split(/\|/, $line);
		if(exists $contigs{$lines[0]}){
			print TGT1 ">$lines[0] $lines[1] ", join(" ", @{$contigs{$lines[0]}}), "\n";
			print TGT1 $seq_obj->seq, "\n";
			#$seqout1->write_seq($seq_obj);
		}elsif(exists $incomplete{$lines[0]}){
			my ($proseq, $plen) = &Find_ORF($seq_obj);
			my $gene = ${$incomplete{$lines[0]}}[0];
			my $qlen = ${$incomplete{$lines[0]}}[2];
			if($mode eq "abs"){
				if($plen > ($qlen-$cutoff)){
					## consider as full assembled contig
					my @alt_contigs = @{$full{$gene}};
					foreach my $alt_contig (@alt_contigs){
						${$incomplete{$alt_contig}}[3] = 1; ## mark all alt contigs to be fully assembled
					}
					print TGT1 ">$lines[0] $lines[1] ", join(" ", @{$incomplete{$lines[0]}}), "\n";
					print TGT1 $seq_obj->seq, "\n";
				}
			}elsif($mode eq "pct"){
				my $percent = $plen/$qlen;
				if($percent > 1-$cutoff){
					## consider as full assembled contig
					my @alt_contigs = @{$full{$gene}};
					foreach my $alt_contig (@alt_contigs){
						${$incomplete{$alt_contig}}[3] = 1; ## mark all alt contigs to be fully assembled
					}
					print TGT1 ">$lines[0] $lines[1] ", join(" ", @{$incomplete{$lines[0]}}), "\n";
					print TGT1 $seq_obj->seq, "\n";
				}
			}else{
				print ERR "ERROR: Please specify correct mode: abs or pct\n";
				die;
			} ## mode if else
		}else{	## the contig doesn't map any gene
			print TGT2 ">$lines[0] $lines[1]", "\n";
			print TGT2 $seq_obj->seq, "\n";
		}## if exists contig has gene mapped
		
	} ## while reading next seq record				
	
	## re-read the assembled contig file, put partially assembled contig into next run
	my $seqio_obj2 = Bio::SeqIO->new(-file => $seqfile, -format => "fasta");
	while(my $seq_obj = $seqio_obj->next_seq){
		my $line = $seq_obj->id;
		my @lines = split(/\|/, $line);
		if(exists $incomplete{$lines[0]}){
			if(${$incomplete{$lines[0]}}[3] == 0){
				print TGT2 ">$lines[0] $lines[1] ", join(" ", @{$incomplete{$lines[0]}}), "\n";
				print TGT2 $seq_obj->seq, "\n";
			}
		}
	}
	
}

close SRC;
close REF;
close TRI;
close TGT1;
close TGT2;
close ERR;
close OUT;

system("grep -E 'ERROR|Error|error' 00.script/shell.script/transfer.saturate.seq.e > 00.script/shell.script/summary.error.log");
system("echo 'success' > 00.script/shell.script/transfer.saturate.seq.log");

system("echo 'Finished 10.transfer.saturate.seq.pl!' >> job.monitor.txt");

########################
sub Find_ORF{
	my $seq_obj = shift @_;
	my @pro_seqs = ();

	## get all six posssible translation
	push @pro_seqs, $seq_obj->translate(-frame => 0)->seq;
	push @pro_seqs, $seq_obj->translate(-frame => 1)->seq;
	push @pro_seqs, $seq_obj->translate(-frame => 2)->seq;
	push @pro_seqs, $seq_obj->revcom->translate(-frame => 0)->seq;
	push @pro_seqs, $seq_obj->revcom->translate(-frame => 1)->seq;
	push @pro_seqs, $seq_obj->revcom->translate(-frame => 2)->seq;

	## find the longest ORF between two stop codons
	my $protein = 0;
	my $prolen = 0;
	foreach my $pro_seq (@pro_seqs){
		my @orfs = split(/\*/, $pro_seq);
		foreach my $orf (@orfs){
			if(length($orf) > $prolen){
				$protein = $orf;
				$prolen = length($orf);
			}
		}
	}

	## chop AAs before the start codon
	$protein =~ s/^[A-KN-Z]+M/M/;
	## add stop codon to the end
	$protein = $protein."*";
	## update protein length
	$prolen = length($protein);
	
	## return value
	return ($protein, $prolen);
}


