#!/usr/bin/perl -w
# run the script:
# time perl 00.script/c4.get.full.genes.pl 01.data/04.GeneOfInterest/GeneID.v2.txt 01.data/05.splitGenes/03.Full.Length 05.Trinity/full.length.contigs.fasta 08.GlobLocal.Comparison

use strict;

my $reffile = shift @ARGV;			# gene information file with gene length, expression, etc. 
my $locfolder = shift @ARGV;		# local assembly folder containing sub-group fully assembled contigs
my $globfile = shift @ARGV;			# trinity assembly file containing assembled contigs
my $outfolder = shift @ARGV;		# output folder for comparison
my $bothout = "$outfolder/both.full.gene.txt";			# output file with common genes
my $locout = "$outfolder/local.uniq.gene.txt";			# output file for unique local genes
my $globout = "$outfolder/global.uniq.gene.txt";			# output file for unique global genes

# build output dir
system("mkdir -p $outfolder");

# get run information
opendir(SRC, $locfolder);
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
# store local fully assembled contigs
my %hash = ();  
open(SRC, "$locfolder/count1");
foreach my $line (<SRC>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	
	my @index = (0, 2..$#lines);	## get other info except gene ID
	my @info = @lines[@index];
	my $run = pop @info;
	
	if(not exists $hash{$lines[1]}){  ## first time to see full length gene
		$hash{$lines[1]} = [$run, \@info];
	}else{							## not the first time
		my @temp = @{$hash{$lines[1]}};
		my $oldrun = $temp[0];
		
		if($oldrun == $run){		## within the same run, contigs correspond to isoform
			push @{$hash{$lines[1]}}, \@info;
		}elsif($run > $oldrun){		## further runs recover the same gene, replace old info
			$hash{$lines[1]} = [$run, \@info];
		}else{
			die "Something unexpected happened.\n";
		}
	}
}
close SRC;

# store global fully assembled contigs
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

## write to output file
open(TGT, ">$bothout");
open(LOT, ">$locout");
open(GOT, ">$globout");

print TGT "Gene\tGlobal\tLocal\tPro_len\tTrans_Len\tGC\tExpr\n";
print LOT "Gene\tLocal\tPro_len\tTrans_Len\tGC\tExpr\n";
print GOT "Gene\tGlobal\tPro_len\tTrans_Len\tGC\tExpr\n";

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
	$gc = sprintf('%.2f', $gc);

	if(not exists $hash{$key}){
		print GOT "$key\t", join(",", @globinfo), "\t$len\t$tranlen\t$gc\t$expr\t", "\n";
		next;
	}
	
	my @local = @{$hash{$key}};
	my $run = shift @local;
	my @locinfo = ();
	foreach my $r (@local){
		my @info = @$r;
		push @locinfo, $info[0].":".$group.":".$run;
	}
	
	print TGT "$key\t", join(",", @globinfo), "\t", join(",", @locinfo), "\t$len\t$tranlen\t$gc\t$expr\t", "\n";
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
	my $run = shift @local;
	my @locinfo = ();
	foreach my $r (@local){
		my @info = @$r;
		push @locinfo, $info[0].":".$group.":".$run;
	}
	$gc = sprintf('%.2f', $gc);
	
	print LOT "$key\t", join(",", @locinfo), "\t$len\t$tranlen\t$gc\t$expr\t", "\n";
}

close TGT;
close LOT;
close GOT;