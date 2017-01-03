#!/usr/bin/perl -w
# run the script: time 00.script/c13.compare.Trinity.Local.pl 08.full.length/ptr.full.contigs.fasta 08.full.length/Final.v2.ptr.blastx.out 05.Trinity/ptr.full.contigs.fasta 05.Trinity/Trinity.new.ptr.blastx.out 12.Local.Trinity.comparison 01.data/00.PriorData/ptr.proteome.fa

use strict;
use Bio::SeqIO;

my $local_seq_file = shift @ARGV;
my $local_blast_file = shift @ARGV;
my $trinity_seq_file = shift @ARGV;
my $trinity_blast_file = shift @ARGV;
my $outfolder = shift @ARGV;
my $dbfile = shift @ARGV;

my $local_seq_obj = Bio::SeqIO->new(-file => $local_seq_file, -format => 'fasta');
my $trinity_seq_obj = Bio::SeqIO->new(-file => $trinity_seq_file, -format => 'fasta');

open(TGT1, ">$outfolder/Local_uniq.txt");
open(TGT2, ">$outfolder/Trinity_uniq.txt");
open(TGT3, ">$outfolder/Trinity_Local_intersect.txt");
system("mkdir -p $outfolder");

my %db_length = ();
print "Reading database sequence file $dbfile...\n";
my $db_seq_obj = Bio::SeqIO->new(-file => $dbfile, -format => 'fasta');
while(my $seqio = $db_seq_obj->next_seq){
	my $id = $seqio->id;
	my $seq = $seqio->seq;
	my $len = length($seq);
	
	if(not exists $db_length{$id}){
		$db_length{$id} = $len;
	}
}

my %local_gene = ();
my %local_blast_store = ();
my %trinity_gene = ();
my %trinity_blast_store = ();
#=cut
print "Reading local assembly sequence file...\n";
while(my $seqio = $local_seq_obj->next_seq){
	my $id = $seqio->id;
	my @desc = split(/\s+/, $seqio->description);
	my $seq = $seqio->seq;
	
	if(not exists $local_gene{$id}){
		$local_gene{$id} = { hitname => $desc[1]
		};
	}
}

print "Reading trinity assembly sequence file...\n";
while(my $seqio = $trinity_seq_obj->next_seq){
	my $id = $seqio->id;
	my @desc = split(/\s+/, $seqio->description);
	my $seq = $seqio->seq;
	
	if(not exists $trinity_gene{$id}){
		$trinity_gene{$id} = { hitname => $desc[1]
		};
	}
}

my %local_blast = ();
my %trinity_blast = ();
open(BLT1, "$local_blast_file");
open(BLT2, "$trinity_blast_file");

print "Reading local assembly blast file...\n";
while(my $line = <BLT1>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	my $qname = $lines[0];
	my $hitname = $lines[1];
	my $mismatch = $lines[4];
	my $bits = $lines[11];
	my $full_length = 0;
	
	if(exists $local_gene{$qname}){$full_length = 1;}
	
	if(not exists $local_blast{$qname}){
		$local_blast{$qname} = { record => $line,
								 bits => $bits,
								 status => 'single',
								 mismatch => $mismatch,
								 hitname => $hitname,
								 full_length => $full_length
		};
	}elsif($bits > $local_blast{$qname}->{bits}){
		$local_blast{$qname} = { record => $line,
								 bits => $bits,
								 status => 'multiple',
								 mismatch => $mismatch,
								 hitname => $hitname,
								 full_length => $full_length
								 };
	}
	
	if(exists $local_gene{$qname} and not exists $local_blast_store{$hitname}){
		$local_blast_store{$hitname} = 0;
	}

}

print "Reading trinity assembly blast file...\n";
while(my $line = <BLT2>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	my $qname = $lines[0];
	my $hitname = $lines[1];
	my $mismatch = $lines[4];
	my $bits = $lines[11];
	my $full_length = 0;
	
	if(exists $trinity_gene{$qname}){$full_length = 1;}

	if(not exists $trinity_blast{$qname}){
		$trinity_blast{$qname} = { record => $line,
								 bits => $bits,
								 status => 'single',
								 mismatch => $mismatch,
								 hitname => $hitname,
								 full_length => $full_length
		};
	}elsif($bits > $trinity_blast{$qname}->{bits}){
		$trinity_blast{$qname} = { record => $line,
								 bits => $bits,
								 status => 'multiple',
								 mismatch => $mismatch,
								 hitname => $hitname,
								 full_length => $full_length
								 };
	}
	if(exists $trinity_gene{$qname} and not exists $trinity_blast_store{$hitname}){
		$trinity_blast_store{$hitname} = 0;
	}

}

print "Finding local assembly unique transcripts...\n";
foreach my $contig_id (sort keys %local_gene){
	if(not exists $local_blast{$contig_id}){next;}
	my $hitname = $local_blast{$contig_id}->{hitname};
	#print "$contig_id\t$hitname\n";
	if(not exists $trinity_blast_store{$hitname}){
		#print TGT1 "Local_uniq\t$contig_id\t$hitname\n";
		#if($local_blast{$contig_id}->{mismatch} > 0){next;}
		print TGT1 ">$hitname\t", $db_length{$hitname}, "\n";
		print TGT1 "Local\t", $local_blast{$contig_id}->{record}, "\n";
		
		my $mark = 1;
		foreach my $trinity_id (keys %trinity_blast){
			my $trinity_hitname = $trinity_blast{$trinity_id}->{hitname};
			if($hitname eq $trinity_hitname){
				print TGT1 "Trinity\t", $trinity_blast{$trinity_id}->{record}, "\n";
				$mark = 0;
			}
		}
		if($mark == 1){print TGT1 "NULL\n";}
	}else{
		print TGT3 ">$hitname\t", $db_length{$hitname}, "\n";
		my $mark = 1;
		foreach my $local_id (keys %local_blast){
			my $local_hitname = $local_blast{$local_id}->{hitname};
			if($hitname eq $local_hitname and $local_blast{$local_id}->{full_length}){
				print TGT3 "Local\t", $local_blast{$local_id}->{record}, "\n";
				$mark = 0;
			}
		}
		if($mark == 1){print TGT3 "NULL\n";}
		
		$mark = 1;
		foreach my $trinity_id (keys %trinity_blast){
			my $trinity_hitname = $trinity_blast{$trinity_id}->{hitname};
			if($hitname eq $trinity_hitname and $trinity_blast{$trinity_id}->{full_length}){
				print TGT3 "Trinity\t", $trinity_blast{$trinity_id}->{record}, "\n";
				$mark = 0;
			}
		}
		if($mark == 1){print TGT3 "NULL\n";}
	}
}

print "Finding trinity assembly unique transcripts...\n";
foreach my $contig_id (sort keys %trinity_gene){
	if(not exists $trinity_blast{$contig_id}){next;}
	my $hitname = $trinity_blast{$contig_id}->{hitname};
	if(not exists $local_blast_store{$hitname}){
		#print TGT2 "Trinity_uniq\t$contig_id\t$hitname\n";
		#if($trinity_blast{$contig_id}->{mismatch} > 0){next;}
		print TGT2 ">$hitname\t", $db_length{$hitname}, "\n";
		print TGT2 "Trinity\t", $trinity_blast{$contig_id}->{record}, "\n";
		
		my $mark = 1;
		foreach my $local_id (keys %local_blast){
			my $local_hitname = $local_blast{$local_id}->{hitname};
			if($hitname eq $local_hitname){
				print TGT2 "Local\t", $local_blast{$local_id}->{record}, "\n";
				$mark = 0;
			}
		}
		if($mark == 1){print TGT2 "NULL\n";}
	}
}

print "Done! :)";

#=cut
close BLT1;
close BLT2;
close TGT1;
close TGT2;
close TGT3;