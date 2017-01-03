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

=pod
## check if previous step has succesfully finished
my $reffolder = "01.data/05.splitGenes/01.Protein/run.0";
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[0-9]+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("blastx.back.$chk.sh.o*");
		my @stdout = glob("blastx.back.$chk.sh.e*");
		my @log = glob("00.script/07.blastx.script/run.$run/$chk.done.*");
		if(!($stderr[0] and $stdout[0] and $log[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' blastx.back.$chk.sh.e* >> 00.script/07.blastx.script/run.$run/summary.error.log\n");
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
					#system("rm -f blastx.back.$chk.*");
					#system("rm -f 00.script/07.blastx.script/run.$run/$chk.done.log");
					system("echo 'Resubmitting the job: truncate.header.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/07.blastx.script/run.$run/blastx.back.$chk.sh");
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

=cut

## start running the script
system("mv blastx.back.* 00.script/07.blastx.script/run.$run/");
system("mv blastn.back.* 00.script/07.blastn.script/run.$run/");

system("rm -rf 00.script/10.transfer.script/run.$run");
system("mkdir -p 00.script/10.transfer.script/run.$run");
system("mkdir -p $tgtfolder");

opendir(SRC, $blastfolder) or die "ERROR: Cannot open $blastfolder: $!";
my @subs = sort(grep(/^[0-9]+$/, readdir SRC));
closedir SRC;

open(ERR, ">$errfile") or die "ERROR: Cannot write $errfile: $!";
open(OUT, ">$outfile") or die "ERROR: Cannot write $outfile: $!";
foreach my $sub (@subs){
	print OUT "$sub\n";
	
	my $trinity_obj = Bio::SeqIO->new(-file => "$srcfolder/$sub/Trinity.new.fasta", -format => 'fasta');
	my %query_length = ();
	while(my $seqio = $trinity_obj->next_seq){
		my $id = $seqio->id;
		my $seq = $seqio->seq;
		my $len = length($seq);
		if(not exists $query_length{$id}){
			$query_length{$id} = $len;
		}
	}
	
	open(BLS, "$blastfolder/$sub/$sub.contigs.blast.out") or die "ERROR: Cannot open $blastfolder/$sub/$sub.contigs.blast.out: $!";
	my %valid = ();
	foreach my $line (<BLS>){
		chomp $line;
		my @lines = split(/\t/, $line);
		my $qname = $lines[0];
		my $qlen = $query_length{$qname};
		my $hitname = $lines[1];
		my $qstart = $lines[6];
		my $qend = $lines[7];
		my $hitstart = $lines[8];
		my $hitend = $lines[9];
		my $evalue = $lines[10];
		my $bits = $lines[11];
		my $frame = 0;

		## find open reading frame
		if($qstart < $qend and $hitstart < $hitend){
			if($qstart % 3 == 1){
				$frame = 1;
			}elsif($qstart % 3 == 2){
				$frame = 2;
			}else{
				$frame = 3;
			}
		}elsif($qstart > $qend and $hitstart < $hitend){
			if(($qlen - $qstart) % 3 == 0){
				$frame = -1;
			}elsif(($qlen - $qstart) % 3 == 1){
				$frame = -2;
			}else{
				$frame = -3;
			}
		}
		
		if(not exists $valid{$qname} or ($bits > $valid{$qname}->{bits})){
            $valid{$qname} = {bits => $bits, frame => $frame};
        }
		
	}
	close(BLS);
	
	$trinity_obj = Bio::SeqIO->new(-file => "$srcfolder/$sub/Trinity.new.fasta", -format => 'fasta');
	system("mkdir -p $tgtfolder/$sub");
	open(TGT1, ">$tgtfolder/$sub/$sub.fasta") or die "ERROR: Cannot write $tgtfolder/$sub/$sub.fasta: $!";
	open(TGT2, ">$tgtfolder/$sub/$sub.prot.fasta") or die "ERROR: Cannot write $tgtfolder/$sub/$sub.prot.fasta: $!";
	
	while(my $seqio = $trinity_obj->next_seq){
		my $id = $seqio->id;
		my @desc = split(/\s+/, $seqio->description);
		my $seq = $seqio->seq;
		
		if(exists $valid{$id}){
			print TGT1 ">$id ", $desc[0], "\n";
			print TGT1 "$seq\n";
			
			my ($fullseq, $proseq, $qlen) = &Find_ORF($seqio, $valid{$id}->{frame});
			print TGT2 ">$id $qlen\n";
			print TGT2 "$fullseq\n";
		}
	}
}

close TGT1;
close TGT2;
close ERR;
close OUT;

system("grep -E 'ERROR|Error|error' 00.script/10.transfer.script/run.$run/transfer.saturate.seq.e > 00.script/10.transfer.script/run.$run/summary.error.log");
system("echo 'success' > 00.script/10.transfer.script/run.$run/transfer.saturate.seq.log");

system("echo 'Finished 10.transfer.saturate.seq.pl!' >> job.monitor.txt");

########################
sub Find_ORF{
	my $seq_obj = shift @_;
	my $frame = shift @_;
	my $pro_seq = 0;
	my $protein = 0;

	## get all six posssible translation
	if($frame == 1){$pro_seq = $seq_obj->translate(-frame => 0)->seq;}
	elsif($frame == 2){$pro_seq = $seq_obj->translate(-frame => 1)->seq;}
	elsif($frame == 3){$pro_seq = $seq_obj->translate(-frame => 2)->seq;}
	elsif($frame == -1){$pro_seq = $seq_obj->revcom->translate(-frame => 0)->seq;}
	elsif($frame == -2){$pro_seq = $seq_obj->revcom->translate(-frame => 1)->seq;}
	elsif($frame == -3){$pro_seq = $seq_obj->revcom->translate(-frame => 2)->seq;}

	## find the longest ORF between two stop codons
	my $prolen = 0;
	my $stop = 0;
	my @orfs = split(/\*/, $pro_seq);
	if(scalar @orfs > 1){$stop = 1;}
	foreach my $orf (@orfs){
		if(length($orf) > $prolen){
			$protein = $orf;
			$prolen = length($orf);
		}
	}

	## chop AAs before the start codon
	$protein =~ s/^[A-L|N-Z]+M/M/;
	## add stop codon to the end
	if($stop){$protein = $protein."*";}
	## update protein length
	$prolen = length($protein);
	
	## return value
	return ($pro_seq, $protein, $prolen);
}
