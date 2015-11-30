#!/usr/bin/perl -w
# run the script: time perl 00.script/11.summarize.run.pl 01.data/04.GeneOfInterest/GeneID.v2.txt 11.stable.run/all.summary.run.txt 07.map.back/03.bowtie.nucl 9

use strict;

my $genefile = shift @ARGV;
my $outfile = shift @ARGV;
my $blastfolder = shift @ARGV;
my $maxrun = shift @ARGV;

my @path = split(/\//, $outfile);
pop @path;
my $path = join("/", @path);
system("mkdir -p $path");

open(GEN, $genefile) or die "ERROR: Cannot open $genefile: $!";
open(TGT, ">$outfile") or die "ERROR: Cannot open $outfile: $!";

my %hash = ();
<GEN>;
foreach my $line (<GEN>){
	chomp $line;
	my @lines = split(/\t/, $line);
	my ($gene, $pro_len, $trans_len, $gc, $gene_len, $group, $fpkm) = @lines[0..6];
	if(not exists $hash{$group}{$gene}){
		$hash{$group}{$gene}{"info"} = [$fpkm, $pro_len, $trans_len, $gene_len, $gc];
	}
}

foreach my $run (0..$maxrun){
	print "This is run $run\n";
	foreach my $key (sort(keys %hash)){
		open(SRC, "$blastfolder/run.$run/$key/$key.contigs.blast.out") or die "ERROR: Cannot open $blastfolder/run.$run/$key/$key.contigs.blast.out: $!";
		
		foreach my $line (<SRC>){
			chomp $line;
			my @lines = split(/\t/, $line);
			my ($contig_all, $gene, $identity, $qstart, $qend, $sstart, $send) = (@lines[0..2], @lines[6..9]);
			my @contig = split(/\|/, $contig_all);
			my $contig = $contig[0];
			if(exists $hash{$key}{$gene}){
				if(not exists $hash{$key}{$gene}{"run.$run"}){
					$hash{$key}{$gene}{"run.$run"} = [$qstart, $qend, $sstart, $send, $identity];
				}else{ # redundant hit
					if((abs($sstart-$send) > 0.8*abs(${$hash{$key}{$gene}{"run.$run"}}[2]-${$hash{$key}{$gene}{"run.$run"}}[3])
						and $identity > 0.8*${$hash{$key}{$gene}{"run.$run"}}[4])){
						$hash{$key}{$gene}{"run.$run"} = [$qstart, $qend, $sstart, $send, $identity];
					}
				} # if exist run
			} # if exist gene
		} # foreach blast line
		
		close SRC;
	} # foreach group
} # foreach run

print TGT "Gene_ID\tFPKM\tPro_len\tTrans_len\tGene_len\tFinal_GC\tGroup_ID\tStable_run\tFull_len\tHit\tIdentity\t", join("\t",(0..$maxrun)),"\n";
foreach my $key (sort(keys %hash)){
	foreach my $gene (sort(keys %{$hash{$key}})){
		my $stable = 0;
		my $full = 0;
		my $get = 0;
		my @qlen = ();
		my @slen = ();
		my @frac = ();
		my $identity = 0;
		my $original = ${$hash{$key}{$gene}{"info"}}[1];
		foreach my $run (0..$maxrun){
			my $qlen = 0;
			my $slen = 0;
			if($hash{$key}{$gene}{"run.$run"}){
				$qlen = abs(${$hash{$key}{$gene}{"run.$run"}}[0] - ${$hash{$key}{$gene}{"run.$run"}}[1])+1;
				$slen = abs(${$hash{$key}{$gene}{"run.$run"}}[2] - ${$hash{$key}{$gene}{"run.$run"}}[3])+1;
			}
			push @qlen, $qlen;
			push @slen, $slen;
			push @frac, $slen/$original;
		}
		my $postqlen = pop @qlen;
		my $postslen = pop @slen;
		if($postslen == 0){
			$get = "absent";
		}else{
			$get = "present";
		}
		if($postslen == $original){
			$full = "T";
		}else{
			$full = "F";
		}
		foreach my $run (reverse 1..$maxrun){
			my $qlen=pop @qlen;
			my $slen=pop @slen;
			$stable = $run;
			if($qlen != $postqlen or $slen != $postslen){
				last;
			}
			$postqlen = $qlen;
			$postslen = $slen;
		}
		if($hash{$key}{$gene}{"run.$stable"}){
			$identity = ${$hash{$key}{$gene}{"run.$stable"}}[4];
		}else{
			$identity = "NA";
		}
		print TGT "$gene\t", join("\t", @{$hash{$key}{$gene}{"info"}}), "\t$key\t$stable\t$full\t$get\t$identity\t", join("\t", @frac), "\n";
	}
}

close GEN;
close TGT;