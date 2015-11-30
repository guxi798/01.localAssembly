#!/usr/bin/perl -w
# run the script: time perl 00.script/10.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.0 07.map.back/03.bowtie.nucl/run.0 01.data/05.splitGenes/02.Transcript/run.1 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 

use strict;
use Bio::SeqIO;

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
my ($run) = $srcfolder =~ /\/(run\.[0-9]+)/;

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
	if($line eq ""){next;}
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}){
		$hash{$lines[0]} = $lines[1];		# key: gene ID, value: protein len
	}
}

## need contig length for inferring frame position
open(LEN, "$srcfolder/Trinity.new.fasta");
my %len = ();
foreach my $line (<LEN>){
	if(!($line =~ /^>/)){next;}
	chomp $line;
	$line =~ s/^>//;
	my @lines = split(/\s+/, $line);
	$lines[1] =~ s/len=//;
	if(not exists $len{$lines[0]}){
		$len{$lines[0]} = $lines[1];
	}
}
close LEN;

my $frame = 0;
open(BLS, "$blastfolder/blast.out") or die "ERROR: Cannot open $blastfolder/blast.out: $!";
my %full = ();
foreach my $line (<BLS>){		## loop over blast record
	chomp $line;
	my @lines = split(/\t/, $line);
	my ($contig_info, $gene) = @lines[0..1];
	my @genes = split(/\./, $gene);
	pop @genes;
	$gene = join(".", @genes);
	my @contig_info = split(/\|/, $contig_info);
	my $contig = $contig_info[0];
    if($lines[6]<$lines[7] and $lines[8]<$lines[9]){
		if($lines[6] % 3 == 1){
			$frame = 1;
		}elsif($lines[6] % 3 == 2){
			$frame = 2;
		}else{
 			$frame = 3;
		}
	}elsif($lines[6]>$lines[7] and $lines[8]<$lines[9]){
		my $contig_len = $len{$contig};
		if(($contig_len - $lines[6]) % 3 == 0){
 			$frame = -1;
		}elsif(($contig_len - $lines[6]) % 3 == 1){
			$frame = -2;
		}else{
			$frame = -3;
		}
	}
	my $left = ($lines[8], $lines[9])[$lines[8] > $lines[9]];		## subject start
	my $right = ($lines[8], $lines[9])[$lines[8] < $lines[9]];		## subject end
	my $mark = 0;
	
	if(not exists $full{$gene}){
		$full{$gene} = ();
		if($left == 1 and $right == $hash{$gene}){
			$mark = 1;
		}
		push @{$full{$gene}}, [$mark, abs($right-$left)+1, $hash{$gene}, $contig, $frame];
	}else{
		if($left == 1 and $right == $hash{$gene}){
			$mark = 1;
		}
		push @{$full{$gene}}, [$mark, abs($right-$left)+1, $hash{$gene}, $contig, $frame];			## there are multiple contigs mapped to the same gene
	}
}
close(BLS);

my %contigs = ();
my %incomplete = ();

foreach my $key (sort(keys %full)){
	foreach my $contig_record (@{$full{$key}}){
		my $mark = shift @$contig_record; 
		my $qlen = shift @$contig_record;	## query length
		my $slen = shift @$contig_record;	## subject length
		my $contig = shift @$contig_record;
		my $frame = shift @$contig_record;

		## the contig covers full length protein
		if($mark == 1){			
			## remove all contigs under the same gene, not give chance to alt contig to grow 
			if(not exists $contigs{$contig}){
				$contigs{$contig} = [$key, $qlen, $slen, $mark, $frame];
			}elsif($key eq ${$contigs{$contig}}[0]){
				print ERR "Warnings: $contig is mapped to the same gene with multiple hit.\n";
			}else{
				print ERR "ERROR: something unexpected happened for $contig\n";
			}
		}else{	## the contig covers only partial length protein
			if(not exists $incomplete{$contig}){
				$incomplete{$contig} = [$key, $qlen, $slen, $mark, $frame];
			}elsif($key eq ${$incomplete{$contig}}[0]){
				print ERR "Warnings: $contig is mapped to the same gene with multiple hit.\n";
			}else{
				print ERR "ERROR: something unexpected happened for $contig\n";
			}
		}
	## the contig covers full length protein
	}
}

#open(TRI, "$srcfolder/$sub/Trinity.new.fasta") or die "ERROR: Cannot open $srcfolder/$sub/Trinity.new.fasta: $!";
system("mkdir -p $tgtfolder");
open(TGT2, ">$tgtfolder/incomplete.fasta") or die "ERROR: Cannot write $tgtfolder/incomplete.fasta: $!";
my $seqfile = "$srcfolder/Trinity.new.fasta";

my $seqio_obj = Bio::SeqIO->new(-file => $seqfile, -format => "fasta");

while(my $seq_obj = $seqio_obj->next_seq){
	my $line = $seq_obj->id;
	#print "$line\n";
	my @lines = split(/\s+/, $line);
	if(exists $contigs{$lines[0]}){
		print TGT1 ">$lines[0] ", join(" ", @{$contigs{$lines[0]}}), "\n";
		print TGT1 $seq_obj->seq, "\n";
	}elsif(exists $incomplete{$lines[0]}){
		my $gene = ${$incomplete{$lines[0]}}[0];
		my $qlen = ${$incomplete{$lines[0]}}[2];
		my $mark = ${$incomplete{$lines[0]}}[3];
		my $frame = ${$incomplete{$lines[0]}}[4];
		my ($proseq, $plen) = &Find_ORF($seq_obj, $frame);
		if($mode eq "abs"){
			if($plen > ($qlen-$cutoff)){
				## consider as full assembled contig
				#my @alt_contigs = @{$full{$gene}};
				#foreach my $alt_contig (@alt_contigs){
				#	${$incomplete{$alt_contig}}[3] = 1; ## mark all alt contigs to be fully assembled
				#}
				print TGT1 ">$lines[0] $gene $plen $qlen $mark $frame \n";
				print TGT1 $seq_obj->seq, "\n";
			}else{
				print TGT2 ">$lines[0] $gene $plen $qlen $mark $frame \n";
				print TGT2 $seq_obj->seq, "\n";
			}
		}elsif($mode eq "pct"){
			my $percent = $plen/$qlen;
			if($percent > 1-$cutoff){
				## consider as full assembled contig
				#my @alt_contigs = @{$full{$gene}};
				#foreach my $alt_contig (@alt_contigs){
				#	${$incomplete{$alt_contig}}[3] = 1; ## mark all alt contigs to be fully assembled
				#}
				print TGT1 ">$lines[0] $gene $plen $qlen $mark $frame \n";
				print TGT1 $seq_obj->seq, "\n";
			}else{
				print TGT2 ">$lines[0] $gene $plen $qlen $mark $frame \n";
				print TGT2 $seq_obj->seq, "\n";
			}
		}else{
			print ERR "ERROR: Please specify correct mode: abs or pct\n";
			die;
		} ## mode if else
	}else{	## the contig doesn't map any gene
		next;
	}## if exists contig has gene mapped
	
} ## while reading next seq record				
	

close REF;
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
	return ($protein, $prolen);
}



