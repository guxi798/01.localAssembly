#!/usr/bin/perl -w
# run the script:
# time perl 00.script/c4.get.full.genes.pl ../03.Protein.Subset/01.data/05.splitGenes/03.Full.Length 01.data/04.GeneOfInterest/GeneID.v2.txt 05.Full.Length/full.length.contigs.fasta 05.Full.Length/both.full.gene.txt 05.Full.Length/local.uniq.gene.txt 05.Full.Length/global.uniq.gene.txt

use strict;

my $srcfolder = shift @ARGV;
my $reffile = shift @ARGV;
my $globfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $locout = shift @ARGV;
my $globout = shift @ARGV;

opendir(SRC, $srcfolder);

my @runs = sort grep(/run\.([0-9]+)/, readdir SRC);
map{$_ =~ s/run\.//} @runs;
@runs = sort {$a <=> $b} @runs;

closedir SRC;

## get protein length and group from reference file
## protein length used for identity calculation
## group ID used for retrieve local assembly
open(REF, $reffile);
my %ref = ();
foreach my $line (<REF>){
	chomp $line;
	my @lines = split(/\t/, $line);
	if(not exists $ref{$lines[0]}){
		$ref{$lines[0]} = [$lines[1], $lines[5], $lines[2], $lines[3], $lines[6]];
	}
}
close REF;

## store full length record
my %hash = ();  
foreach my $run (@runs){
	open(SRC, "$srcfolder/run.$run/full.length.contigs.fasta");
	foreach my $line (<SRC>){
		chomp $line;
		if(!($line =~ /^>/)){next;}
		$line =~ s/^>//;
		my @lines = split(/\s+/, $line);
		
		my @index = (0, 2..$#lines);	## get other info except gene ID
		my @info = @lines[@index];
		push @info, $run;		## store run info
		
		if(not exists $hash{$lines[1]}){  ## first time to see full length gene
			$hash{$lines[1]} = [\@info];
		}else{							## not the first ime
			my @record = @{$hash{$lines[1]}};
			my @oldinfo = @{$record[0]};
			my $oldrun = pop @oldinfo;
			
			if($oldrun == $run){		## within the same run, contigs correspond to isoform
				push @{$hash{$lines[1]}}, \@info;
			}elsif($oldrun < $run){		## further runs recover the same gene, replace old info
				$hash{$lines[1]} = [\@info];
			}else{
				die "Something unexpected happened.\n";
			}
		}
	}
	close SRC;
}

open(GBL, $globfile);
my %hash2 = ();
foreach my $line (<GBL>){
	chomp $line;
	if($line !~ /^>/){next;}
	$line =~ s/^>//;
	my @lines = split(/\s+/, $line);
	#print join("\t", @lines), "\n";
	
	my @index = (0, 2..$#lines);	## get other info except gene ID
	my @info = @lines[@index];
	
	if(not exists $hash2{$lines[1]}){
		$hash2{$lines[1]} = [\@info];
	}else{
		push @{$hash2{$lines[1]}}, \@info;
	}
}
close GBL;

open(TGT, ">$tgtfile");
open(LOT, ">$locout");
open(GOT, ">$globout");
foreach my $key (sort keys %hash2){
	my @global = @{$hash2{$key}};
	my @globinfo = ();
	foreach my $r (@global){
		my @info = @$r;
		push @globinfo, $info[0];
	}

	my $len = 0;
	my $group = 0;
	my $tranlen = 0;
	my $gc = 0;
	my $expr = 0;
	if(exists $ref{$key}){
		my @temp = @{$ref{$key}};
		$len = shift @temp;
		$group = shift @temp;
		$tranlen = shift @temp;
		$gc = shift @temp;
		$expr = shift @temp;
	}else{
		die "Error: not found $key in reference file $reffile\n";
	}

	if(not exists $hash{$key}){
		print GOT "$key\t$len\t$tranlen\t$gc\t$expr\t", join(",", @globinfo), "\n";
		next;
	}
	
	my @local = @{$hash{$key}};
	my @locinfo = ();
	foreach my $r (@local){
		my @info = @$r;
		my $run = pop @info;
		push @locinfo, $info[0].":".$group.":".$run;
	}
	
	print TGT "$key\t$len\t$tranlen\t$gc\t$expr\t", join(",", @globinfo), "\t", join(",", @locinfo), "\n";
}

foreach my $key (sort keys %hash){
	if(exists $hash2{$key}){next;}
	my $len = 0;
	my $group = 0;
	my $tranlen = 0;
	my $gc = 0;
	my $expr = 0;
	if(exists $ref{$key}){
		my @temp = @{$ref{$key}};
		$len = shift @temp;
		$group = shift @temp;
		$tranlen = shift @temp;
		$gc = shift @temp;
		$expr = shift @temp;
	}else{
		die "Error: not found $key in reference file $reffile\n";
	}
	
	my @local = @{$hash{$key}};
	my @locinfo = ();
	foreach my $r (@local){
		my @info = @$r;
		my $run = pop @info;
		push @locinfo, $info[0].":".$group.":".$run;
	}
	
	print LOT "$key\t$len\t$tranlen\t$gc\t$expr\t", join(",", @locinfo), "\n";
}

close TGT;
close LOT;
close GOT;