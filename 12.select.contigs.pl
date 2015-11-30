#!/usr/bin/perl -w
# run the script: time perl 00.script/10.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.0/01.trinity 07.map.back/03.bowtie.nucl/run.0 01.data/05.splitGenes/02.Transcript/run.1 01.data/04.GeneOfInterest/GeneID.v1.txt 

use strict;
use Bio::SeqIO;

system("echo 'Running 10.transfer.saturate.seq.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;			## assembled contig folder
my $blastfolder = shift @ARGV;			## blast result folder
my $tgtfolder = shift @ARGV;			## target folder to put contigs for new run
my $reffile = shift @ARGV;				## reference file with gene ID and protein length info
#my $mode = shift @ARGV;					## abs: absolute value; pct: percent value
#my $cutoff = shift @ARGV;				## absolute AA number, or percent
my $errfile = "00.script/shell.script/transfer.saturate.seq.e";
my $outfile = "00.script/shell.script/transfer.saturate.seq.o";
my ($run) = $srcfolder =~ /\/(run\.[0-9]+)\//;

my $identity_cut = 35
my $len_cut = 0.5

## record file: error and output file
open(ERR, ">$errfile") or die "ERROR: Cannot write $errfile: $!";
open(OUT, ">$outfile") or die "ERROR: Cannot write $outfile: $!";

open(TGT1, ">$tgtfolder/contigs.good.fasta") or die "ERROR: Cannot write $tgtfolder/contigs.good.fasta: $!";
open(TGT2, ">$tgtfolder/contigs.poor.fasta") or die "ERROR: Cannot write $tgtfolder/contigs.poor.fasta: $!";

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
	my %good = ();
	my %poor = ();
	foreach my $line (<BLS>){		## loop over blast record
		chomp $line;
		my @lines = split(/\t/, $line);
		my ($contig_info, $gene, $identity, $map_len) = @lines[0..3];
		my @contig_info = split(/\|/, $contig_info);
		my $contig = $contig_info[0];
		
		if($identity < $identity_cut){
			$poor{$contig} = 0;
			next;
		}
		
		if(exists $hash{$gene}){
			my $pro_len = $hash{$gene}
			if($map_len < $pro_len * $len_cut){
				$poor{$contig} = 0;
				next;
			}
			
			if(not exists $good{$gene}){
				$good{$gene}{$contig} = $map_len;
			}else{
				my @contig1 = keys(%{$good{$gene}});
				foreach my $contig1 (@contig1){
					my $map_len1 = $good{$gene}{$contig1};
					my @part1 = split(/_/, $contig1);
					my @part2 = split(/_/, $contig);
					
					if($part1[0] eq $part2[0]){ ## the same gene, diff isoform
						if ($map_len > $map_len1 and abs($map_len - $map_len1) > $pro_len*0.01){
							delete $good{$gene}{$contig1};
							$good{$gene}{$contig} = $map_len;
						}elsif($contig le $contig1){
							delete $$good{$gene}{$contig1};
							$good{$gene}{$contig} = $map_len;
						}else{
							$poor{$contig} = 0;
						}
					}else{	## diff group
						$good{$gene}{$contig} = $map_len;
					}
					
				}
				
			}
		}else{
			die "Error: $gene doesn't exist in the reference file\n";
		}
	}
	close(BLS);
	
	my %contigs = ();
	my %incomplete = ();
	
	foreach my $key (sort(keys %good)){
		my @contigs = sort(keys %{$good{$key}});
		foreach my $contigs (@contigs){
			$contigs{$contigs} = $key;
		}
	}
	
	#open(TRI, "$srcfolder/$sub/Trinity.new.fasta") or die "ERROR: Cannot open $srcfolder/$sub/Trinity.new.fasta: $!";
	my $seqfile = "$srcfolder/$sub/Trinity.new.fasta";
	
	my $seqio_obj = Bio::SeqIO->new(-file => $seqfile, -format => "fasta");
	
	while(my $seq_obj = $seqio_obj->next_seq){
		my $line = $seq_obj->id;
		#print "$line\n";
		my @lines = split(/\|/, $line);
		
		if(exists $contigs{$lines[0]}){
			print TGT1 ">$lines[0]#$contigs{$lines[0]}\n";
			print TGT1 $seq_obj->seq, "\n";
		}elsif(exists $poor{$lines[0]}){
			print TGT2 ">$lines[0] $lines[1]", "\n";
			print TGT2 $seq_obj->seq, "\n";
		}## if exists contig has gene mapped
		
	} ## while reading next seq record				
}

close SRC;
close REF;
close TRI;
close TGT1;
close TGT2;
close ERR;
close OUT;

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


