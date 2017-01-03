#!/usr/bin/perl
# run the script: time perl 00.script/c14.retrieveseq.Trinity.Local.pl 12.Local.Trinity.comparison/Trinity_Local_results.txt 13.Local.Trinity.alignment 
# 01.data/06.TargetTranscriptome/transcriptome.v1.fa 01.data/06.TargetTranscriptome/proteome.v1.fa 
# 05.Trinity/Trinity.new.fasta 05.Trinity/Trinity.prot.fasta 01.data/00.PriorData/ptr.proteome.fa

use strict;
use Bio::SeqIO;

my $srcfile = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $local_nucl = shift @ARGV;
my $local_prot = shift @ARGV;
my $trinity_nucl = shift @ARGV;
my $trinity_prot = shift @ARGV;
my $ref_prot = shift @ARGV;

my %local_seq = ();
my %trinity_seq = ();
my %ref_seq = ();

system("mkdir -p $tgtfolder");

my $local_nucl_obj = Bio::SeqIO->new(-file => $local_nucl, -format => 'fasta');
my $local_prot_obj = Bio::SeqIO->new(-file => $local_prot, -format => 'fasta');
my $trinity_nucl_obj = Bio::SeqIO->new(-file => $trinity_nucl, -format => 'fasta');
my $trinity_prot_obj = Bio::SeqIO->new(-file => $trinity_prot, -format => 'fasta');
my $ref_prot_obj = Bio::SeqIO->new(-file => $ref_prot, -format => 'fasta');

while(my $seqio = $local_nucl_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $local_seq{$id}{'nucl'}){
		$local_seq{$id}{'nucl'} = $seqio->seq;
	}
}
while(my $seqio = $local_prot_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $local_seq{$id}{'prot'}){
		$local_seq{$id}{'prot'} = $seqio->seq;
	}
}

while(my $seqio = $trinity_nucl_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $trinity_seq{$id}{'nucl'}){
		$trinity_seq{$id}{'nucl'} = $seqio->seq;
	}
}
while(my $seqio = $trinity_prot_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $trinity_seq{$id}{'prot'}){
		$trinity_seq{$id}{'prot'} = $seqio->seq;
	}
}

while(my $seqio = $ref_prot_obj->next_seq){
	my $id = $seqio->id;
	if(not exists $ref_seq{$id}{'prot'}){
		$ref_seq{$id}{'prot'} = $seqio->seq;
	}
}

my $i = 0;
open(SRC, $srcfile);
while(my $line = <SRC>){
	chomp $line;
	if($line =~ /^>/){
		$line =~ s/>//g;
		if($i){close NUC; close PRO;}
		open(NUC, ">$tgtfolder/$line.nucl.fasta");
		open(PRO, ">$tgtfolder/$line.prot.fasta");
	}else{
		$i ++;
		my @lines = split(/\s+/, $line);
		if(exists $trinity_seq{$lines[0]}){
			print NUC ">Trinity.$lines[0]\n", $trinity_seq{$lines[0]}{'nucl'}, "\n";
		}else{
			print "Warning: Not existing $lines[0] in Trinity\n";
		}
		if(exists $local_seq{$lines[1]}){
			print NUC ">Local.$lines[1]\n", $local_seq{$lines[1]}{'nucl'}, "\n";
		}else{
			print "Warning: Not existing $lines[1] in Local\n";
		}
		
		$line = <SRC>;
		my @lines = split(/\s+/, $line);
		if(exists $trinity_seq{$lines[0]}){
			print PRO ">Trinity.$lines[0]\n", $trinity_seq{$lines[0]}{'prot'}, "\n";
		}else{
			print "Warning: Not existing $lines[0] in Trinity\n";
		}
		if(exists $ref_seq{$lines[1]}){
			print PRO ">$lines[1]\n", $ref_seq{$lines[1]}{'prot'}, "\n";
		}
		my $trinity_gene = $lines[1];
		
		$line = <SRC>;
		my @lines = split(/\s+/, $line);
		if(exists $local_seq{$lines[0]}){
			print PRO ">Local.$lines[0]\n", $local_seq{$lines[0]}{'prot'}, "\n";
		}else{
			print "Warning: Not existing $lines[0] in Local\n";
		}
		if(exists $ref_seq{$lines[1]} and $lines[1] ne $trinity_gene){
			print PRO ">$lines[1]\n", $ref_seq{$lines[1]}{'prot'}, "\n";
		}
		my $trinity_gene = $lines[1];
	}
}

close SRC;
close NUC;
close PRO;

