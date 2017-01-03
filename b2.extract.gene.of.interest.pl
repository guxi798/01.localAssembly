#!/usr/bin/perl -w
# run the script: time perl 00.script/b2.extract.gene.of.interest.pl 
#					01.data/00.PriorData/SUT.txt 
#					07.map.back/03.bowtie.nucl 
#					06.assembly/03.bowtie.nucl 
#					01.data/04.GeneOfInterest/GeneID.v1.txt
#					11.stable.run/SUT.blast.summary.txt
#					11.stable.run/SUT.contig.nuclseq.fasta
#					11.stable.run/SUT.contig.protseq.fasta

use strict;
use Bio::SeqIO;

my $reffile = shift @ARGV;			# file contains genes of interests. e.g. SUT
my $blastfolder = shift @ARGV;		# blast result folder
my $contigfolder = shift @ARGV;	# trinity result folder
my $infofile = shift @ARGV;			# gene info file contains protein length and group ID
my $tgtfile1 = shift @ARGV; 		# overall summary output
my $tgtfile2 = shift @ARGV; 		# include detailed sequences
my $tgtfile3 = shift @ARGV; 		# include detailed sequences
my $trinity_blast = shift @ARGV;
my $trinity_nucl = shift @ARGV;
my $trinity_prot = shift @ARGV;

my @path = split(/\//, $tgtfile1);
pop @path;
my $path = join("/", @path);
system("mkdir -p $path");

open(REF, $reffile) or die "ERROR: Cannot open $reffile: $!";
opendir(SRC, $blastfolder) or die "ERROR: Cannot open $blastfolder: $!";
open(INF, $infofile) or die "ERROR: Cannot open $infofile: $!";

#### based on info file, construct a hash reference
my %hash = ();
foreach my $line (<INF>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[0]}){
		$hash{$lines[0]} = [$lines[1], $lines[5], $lines[6]]; # protein length, group ID, expression
	}
}

#### get how many runs
my @runs = grep(/run\.[0-9]+$/, readdir(SRC));
@runs = map {$_ =~ s/run\.//g; $_} @runs;
@runs = sort {$a <=> $b} @runs;
#my @runs = ("run.9");
#### get blast result records for each gene in each run
system("> $tgtfile1");
foreach my $line (<REF>){	# for each GOI
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $hash{$lines[1]}){
		next;
	}
	my $prolen = ${$hash{$lines[1]}}[0];
	my $group = ${$hash{$lines[1]}}[1];
	#my $expr = ${$hash{$lines[1]}}[2];
	#system("echo \">$lines[1]	$lines[1]	ProLen=$prolen	$group $expr\" >> $tgtfile1");
	system("echo \">$lines[1]	$lines[0]	ProLen=$prolen	$group\" >> $tgtfile1");
	system("echo \"\@Trinity\" >> $tgtfile1");
	system("grep \"$lines[1]\" $trinity_blast >> $tgtfile1");
	foreach my $run (@runs){	# for each run
		system("echo \"\@run.$run\" >> $tgtfile1");
		system("grep \"$lines[1]\" $blastfolder/run.$run/$group/$group.contigs.blast.out >> $tgtfile1");
	}
}

#### get contig name from the output file from previous step
open(SUM, $tgtfile1) or die "ERROR: Cannot open $tgtfile1: $!";
#### open sequence output file
open(TGT1, ">$tgtfile2") or die "ERROR: Cannot write $tgtfile2: $!";
open(TGT2, ">$tgtfile3") or die "ERROR: Cannot write $tgtfile3: $!";

my $group = 0;		# group ID
my $run = 0;		# run ID
my $contig = 0;		# contig ID
my %protein = ();	# hash to store contig seqs
my %transcript = ();
my $prev = 0;
my $prun = 0;

my %trinity_transcript = ();
my %trinity_protein = ();
my $trinity_nucl_obj = Bio::SeqIO->new(-file => "$trinity_nucl", -format => 'fasta');
while(my $seqio = $trinity_nucl_obj->next_seq){
	my $id = $seqio->id;
	my $seq = $seqio->seq;
	if(not exists $trinity_transcript{$id}){
		$trinity_transcript{$id} = $seq;
	}
}
my $trinity_prot_obj = Bio::SeqIO->new(-file => "$trinity_prot", -format => 'fasta');
while(my $seqio = $trinity_prot_obj->next_seq){
	my $id = $seqio->id;
	my $seq = $seqio->seq;
	if(not exists $trinity_protein{$id}){
		$trinity_protein{$id} = $seq;
	}
}

my $mode = 0;
foreach my $line (<SUM>){	# for each blast record, corresponding to each contig
	chomp $line;
	my @lines = split(/\s+/, $line);
	if($line =~ /^>/){		# if it is the info line about GOI
		$group = $lines[3];	# get group ID to open corresponding file
		print TGT1 "$line\n";
		print TGT2 "$line\n";
	}elsif($line =~ /^\@run/){	# if it is the info line about run ID
		$mode = 'local';
		$run = $line;
		$run =~ s/\@run\.//;
		print TGT1 "\@run.$run\n";
		print TGT2 "\@run.$run\n";
		$contig = 0;		# initialize contig structure to store contig seqs
		
		# open contig file, that is the trinity output file, group ID and run ID are needed.
		%protein = ();
		%transcript = ();
		my $temp = $run+1;
		
		my $nucl_obj = Bio::SeqIO->new(-file => "$contigfolder/run.$temp/$group/$group.fasta", -format => 'fasta');
		while(my $seqio = $nucl_obj->next_seq){
			my $id = $seqio->id;
			my $seq = $seqio->seq;
			if(not exists $transcript{$id}){
				$transcript{$id} = $seq;
			}
		}
		my $prot_obj = Bio::SeqIO->new(-file => "$contigfolder/run.$temp/$group/$group.prot.fasta", -format => 'fasta');
		while(my $seqio = $prot_obj->next_seq){
			my $id = $seqio->id;
			my $seq = $seqio->seq;
			if(not exists $protein{$id}){
				$protein{$id} = $seq;
			}
		}
	}elsif($line =~ /^\@Trinity/){
		$mode = 'trinity';
		print TGT1 "\@Trinity\n";
		print TGT2 "\@Trinity\n";
	}else{					# if it is the blast record line
		my @ids = split(/\|/, $lines[0]);	# get contig ID
		my $id = $ids[0];
		if($mode eq 'local'){
			if((exists $transcript{$id}) and !($id eq $prev and $prun eq $run)){
				print TGT1 ">$id\n";
				print TGT1 $transcript{$id}, "\n";
			}
			if((exists $protein{$id}) and !($id eq $prev and $prun eq $run)){
				print TGT2 ">$id\n";
				print TGT2 $protein{$id}, "\n";
			}
		}elsif($mode eq 'trinity'){
			if((exists $trinity_transcript{$id}) and !($id eq $prev)){
				print TGT1 ">$id\n";
				print TGT1 $trinity_transcript{$id}, "\n";
			}
			if((exists $trinity_protein{$id}) and !($id eq $prev)){
				print TGT2 ">$id\n";
				print TGT2 $trinity_protein{$id}, "\n";
			}
		}
		$prev = $id;
		$prun = $run;
	}
}

close REF;
close SRC;
close INF;
close TGT1;
close TGT2;
