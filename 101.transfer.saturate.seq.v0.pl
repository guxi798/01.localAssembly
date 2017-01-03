#!/usr/bin/perl -w
# run the script: time perl 00.script/101.transfer.saturate.seq.pl 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/blastx.out 10.unmapped.reads.trinity 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 

use strict;

use Bio::SeqIO;
use Bio::SearchIO;

## read in parameters required by the script
my $usage = "usage: $0 blast+.outfmt5.xml query.fasta output_folder [output_prefix=NameOfBlastFileHere] [mode=abs|pct] [cutoff=NumericCutoff]\n\n";

my $blastfile = shift @ARGV or die $usage;              ## blast result
my $srcfile = shift @ARGV or die $usage;				## assembled contig file
my $tgtfolder = shift @ARGV or die $usage;              ## target folder
my $output_prefix = shift @ARGV || "$blastfile";
my $mode = shift @ARGV || "pct";					## abs: absolute value; pct: percent value
my $cutoff = shift @ARGV || "0.02";                 ## absolute AA number, or percent

system("mkdir -p $tgtfolder");
my $logfile = "$tgtfolder/detect.full.length.log";

open(TGT1, ">$tgtfolder/$output_prefix.full.contigs.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.full.length.contigs.fasta: $!";
open(TGT2, ">$tgtfolder/$output_prefix.incomplete.contigs.fasta") or die "ERROR: Cannot write $tgtfolder/$output_prefix.incomplete.fasta: $!";
open(TGT3, ">$tgtfolder/$output_prefix.unmapped.contigs.fasta");

## read in blast results
my $blast = new Bio::SearchIO(-format => 'blastxml', -file => $blastfile);

my %full_length = ();
my %incomplete_length = ();

while(my $result = $blast->next_result){
    my @names = split(/\s+/, $result->query_description);
    my $qname = $names[0];
    my $qlen = $result->query_length;
	my $frame = 0;
    
    while(my $hit = $result->next_hit){
        my $hitname = $hit->name;
        my $hitlen = $hit->length;
        
        my $qstart = $hit->start('query');
        my $qend = $hit->end('query');
        my $hitstart = $hit->start('hit');
        my $hitend = $hit->end('hit');
		my $qstr = $hit->strand('query');
		
		my $num_hsps = $hit->num_hsps;
		
		while(my $hsp = $hit->next_hsp){
			my $frame2 = $hsp->query->frame;
			my $alen = $hsp->length('query')/3;
			my $evalue = $hsp->evalue;
			my $bits = $hsp->bits;
			
			## find open reading frame
			if($qstr == 1){
				if($qstart % 3 == 1){
					$frame = 1;
				}elsif($qstart % 3 == 2){
					$frame = 2;
				}else{
					$frame = 3;
				}
			}elsif($qstr == -1){
				if(($qlen - $qend) % 3 == 0){
					$frame = -1;
				}elsif(($qlen - $qend) % 3 == 1){
					$frame = -2;
				}else{
					$frame = -3;
				}
			}
			
			#print "$qname\t$qlen\t$hitname\t$hitlen\t$alen\t$qstart\t$qend\t$frame\t$frame2\t$bits\n";

			my $left = ($hitstart, $hitend)[$hitstart > $hitend];		## subject start
			my $right = ($hitstart, $hitend)[$hitstart < $hitend];		## subject end
			my $mark = 0;
			
			if($left == 1 and $right == $hitlen){
				if(not exists $full_length{$qname} or ($bits > $full_length{$qname}->{bits})){
					$full_length{$qname} = { query_id => $qname,
                                             hit_id => $hitname,
											 query_len => $qlen,
											 hit_len => $hitlen,
                                             query_start => $qstart,
                                             query_end => $qend,
                                             hit_start => $hitstart,
                                             hit_end => $hitend,
                                             evalue => $evalue,
                                             bits => $bits,
											 alignment_len => $alen,
											 frame => $frame
					};
					print "$qname\t$qlen\t$hitname\t$hitlen\t$alen\t$qstart\t$qend\t$frame\n";
					if(exists $full_length{$qname}){
						#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
					}
				}else{
					#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
				}
				
			}else{
				if(not exists $incomplete_length{$qname} or ($bits > $incomplete_length{$qname}->{bits})){
					$incomplete_length{$qname} = { query_id => $qname,
												   hit_id => $hitname,
												   query_len => $qlen,
												   hit_len => $hitlen,
												   query_start => $qstart,
												   query_end => $qend,
												   hit_start => $hitstart,
												   hit_end => $hitend,
												   evalue => $evalue,
												   bits => $bits,
												   alignment_len => $alen,
												   frame => $frame
												   };
					print "$qname\t$qlen\t$hitname\t$hitlen\t$alen\t$qstart\t$qend\t$frame\n";
					if(exists $incomplete_length{$qname}){
						#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
					}
				}else{
					#print "Warnings: $qname is mapped to the same gene with multiple hit.\n";
				}
			}        
    }

		}
       
}

my $seqio_obj = Bio::SeqIO->new(-file => $srcfile, -format => "fasta");

while(my $seq_obj = $seqio_obj->next_seq){
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my @temp = split(/\s+/, $seq_obj->desc);
	my $seqlen = $temp[0];
	
	if(exists $full_length{$id}){
		print TGT1 ">$id len=", $full_length{$id}->{query_len}, " ", 
					$full_length{$id}->{hit_id}, " ",
					$full_length{$id}->{hit_len}, " ", 
					$full_length{$id}->{alignment_len}, " blast_full ", 
					$full_length{$id}->{frame}, "\n";
		print TGT1 "$seq\n";
	}elsif(exists $incomplete_length{$id}){
		my ($proseq, $qlen) = &Find_ORF($seq_obj, $incomplete_length{$id}->{frame});
		my $hitlen = $incomplete_length{$id}->{hit_len};
		
		if(($mode eq "abs" and $qlen > ($hitlen-$cutoff) || ($mode eq "pct" and $qlen/$hitlen > 1-$cutoff))){
				## consider as full assembled contig
				print TGT1 ">$id len=", $incomplete_length{$id}->{query_len}, " ", 
						$incomplete_length{$id}->{hit_id}, " ",
						$incomplete_length{$id}->{hit_len}, " ", 
						$qlen, " infer_full ", 
						$incomplete_length{$id}->{frame}, "\n";
				print TGT1 "$seq\n";
		}elsif($mode eq "abs" or $mode eq "pct"){
				print TGT2 ">$id len=", $incomplete_length{$id}->{query_len}, " ", 
						$incomplete_length{$id}->{hit_id}, " ",
						$incomplete_length{$id}->{hit_len}, " ", 
						$qlen, " infer_full ", 
						$incomplete_length{$id}->{frame}, "\n";
				print TGT2 "$seq\n";
		}else{
			print "ERROR: Please specify correct mode: abs or pct\n";
			die;
		} ## mode if else
	}else{
		print TGT3 ">$id $seqlen\n";
		print TGT3 "$seq\n";
	}
}

close TGT1;
close TGT2;
close TGT3;

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



