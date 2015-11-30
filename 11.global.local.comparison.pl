#!/usr/bin/perl -w
# module load perl/5.20.2
# time perl 00.script/11.global.local.comparison.pl 05.Full.Length/both.full.gene.txt 04.blastn/blast.xml.out ../03.Protein.Subset/07.map.back/02.blastn 06.more.detail/both.hsp.num.txt 06.more.detail/both.identity.txt 

use strict;
use Bio::SearchIO;
use List::Util qw(sum);

my $fullfile = shift @ARGV;
my $globlastfile = shift @ARGV;
my $locblastfolder = shift @ARGV;
my $hspfile = shift @ARGV;
my $tgtfile = shift @ARGV;

## put fully recovered genes into hash table
open(FUL, $fullfile);

my %globfull = ();
my %locfull = ();
my %genefull = ();

foreach my $line(<FUL>){
	chomp $line;
	my @lines = split(/\t/, $line);
	if(not exists $genefull{$lines[0]}){
		$genefull{$lines[0]} = 0;
	}else{
		print "This gene appeared more than once: NEED CHECK $lines[0]!\n";
	}
	
	my @glob = split(/,/, $lines[5]);
	foreach my $glob (@glob){
		if(not exists $globfull{$glob}){
			$globfull{$glob} = 0;
		}
	}
	
	my @loc = split(/,/, $lines[6]);
	foreach my $temp (@loc){
		my @temp = split(/:/, $temp);
		my $loc = shift @temp;
		my $group = shift @temp;
		my $run = shift @temp;
		
		if(not exists $locfull{$run}{$group}{$loc}){
			$locfull{$run}{$group}{$loc} = 0;
		}
	}
}
close FUL;

## get in global blast record
print "global blast.....\n";
my $globin = new Bio::SearchIO(-format 	=> 'blastxml',
								-file	=> $globlastfile);
my %globblast = ();
my %globhsps = ();
my %globcomp = ();
my $globcount = 0;

while(my $result = $globin->next_result){
	my @names = split(/\s+/, $result->query_description);
	my $len = $names[1];
	$len =~ s/len=//;
	#if(not exists $globfull{$names[0]}){next;}
	while(my $hit = $result->next_hit){
		$globcount++;
		if($globcount %5000 == 0){print "$globcount\n";}
		my @hitname = split(/\./, $hit->name);
		my $hitname = join(".", @hitname);
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		my $hspnum = $hit->num_hsps;

		my $identity2 = 0;
		my $mark = 'N';
		while(my $hsp = $hit->next_hsp){
			$identity2 += $hsp->num_identical;
		}
		if($identity2 > $identity){
			$mark = 'Y';
		}
		my @info = ($names[0], $hspnum, $mark, $len);
		if(not exists $globhsps{$hspnum}){
			$globhsps{$hspnum} = 1;
		}else{
			$globhsps{$hspnum} ++;
		}
		
		if(not exists $globcomp{$hitname}){
			$globcomp{$hitname} = [\@info];
		}else{
			push @{$globcomp{$hitname}}, \@info;
		}
		#print $identity, "\t",$identity/$len,"\n";
		if(not exists $globblast{$hitname}){
			$globblast{$hitname}{'I'} = [$identity/$len];
			$globblast{$hitname}{'L'} = [$len];
		}else{
			push @{$globblast{$hitname}{'I'}}, $identity/$len;
			push @{$globblast{$hitname}{'L'}}, $len;
		}
	}
}

print "Total: $globcount\n";
foreach my $key (sort {$a <=> $b} keys %globhsps){
	print "Number of hsps: $key\t$globhsps{$key}\n";
}

## get in local blast record
print "local blast.....\n";
my @runs = (12);
my %locblast = ();
my %lochsps = ();
my %loccomp = ();
my $loccount = 0;

foreach my $run (@runs){
	print "\t$locblastfolder/run.$run\n";
	opendir(SUB, "$locblastfolder/run.$run");
	my @groups = sort {$a <=> $b} grep(/[0-9]+/, readdir SUB);

	foreach my $group (@groups){
		print "run.$run\t$group\n";
		my $locin = new Bio::SearchIO(-format => 'blastxml', 
				-file => "$locblastfolder/run.$run/$group/$group.contigs.blast.xml.out");
		my $fullrun = $run + 1;
		while(my $result = $locin->next_result){
			my @names = split(/\s+|\|/, $result->query_description);
			#if(not exists $locfull{$fullrun}{$group}{$names[0]}){next;}
			my $len = $names[1];
			$len =~ s/len=//;
			while(my $hit = $result->next_hit){
				$loccount ++;
				if($loccount % 5000 == 0){print "$loccount\n";}
				my $hitname = $hit->name;
				my $identity = scalar($hit->seq_inds('query', 'identical'));
				my $hspnum = $hit->num_hsps;
				
				my $identity2 = 0;
				my $mark = 'N';
				while(my $hsp = $hit->next_hsp){
					$identity2 += $hsp->num_identical;
				}
				if($identity2 > $identity){
					$mark = 'Y';
				}
				my @info = ($names[0], $hspnum, $mark, $len);
				
				#print $hit->ambiguous_aln() if($hit->ambiguous_aln() ne '-');
				if(not exists $lochsps{$hspnum}){
					$lochsps{$hspnum} = 1;
				}else{
					$lochsps{$hspnum} ++;
				}
				
				if(not exists $loccomp{$hitname}){
					$loccomp{$hitname} = [\@info];
				}else{
					push @{$loccomp{$hitname}}, \@info;
				}
				if(not exists $locblast{$hitname}{$names[0]}){
					$locblast{$hitname}{$names[0]}{'I'} = [$identity/$len];
					$locblast{$hitname}{$names[0]}{'L'} = $len;
				}else{
					push @{$locblast{$hitname}{$names[0]}{'I'}}, $identity/$len;
					push @{$locblast{$hitname}{$names[0]}{'L'}}, $len;
					print "Warning: $hitname + $names[0] has occurred before\n";
				}
			}
		}
	}
	print "\n";
	
	closedir SUB;
}

print "Total: $loccount\n";
foreach my $key (sort {$a <=> $b} keys %lochsps){
	print "Number of hsps: $key\t$lochsps{$key}\n";
}

## write to info output file
open(HSP, ">$hspfile");
foreach my $key (sort keys %globcomp){
	if(exists $loccomp{$key}){
		print HSP "$key\t";
		foreach my $info (@{$globcomp{$key}}){
			my @info = @$info;
			print HSP join(",", @info), ";";
		}
		print HSP "\t";
		foreach my $info (@{$loccomp{$key}}){
			my @info = @$info;
			print HSP join(",", @info), ";";
		}
		print HSP "\n";
	}
}
close HSP;

## re-read fully recovered genes 
open(TGT, ">$tgtfile");
foreach my $key (sort(keys %globblast)){
	if(not exists $locblast{$key}){next;}
	my $globident = 0;
	my $globlen = 0;
	my $locident = 0;
	my $loclen = 0;
	
	my @glob = @{$globblast{$key}{'I'}};
	my @len = @{$globblast{$key}{'L'}};
	$globident = &mean(@glob);
	$globlen = &mean(@len);
	
	my @temp = ();
	@len = ();
	foreach my $contig (sort(keys %{$locblast{$key}})){
		my @loc = @{$locblast{$key}{$contig}{'I'}};
		push @temp, &mean(@loc);
		push @len, $locblast{$key}{$contig}{'L'};
	}
	$locident = &mean(@temp);
	$loclen = &mean(@len);
	
	print TGT "$key\t$globlen\t$globident\t$loclen\t$locident\n";
}
close TGT;

############################
sub mean {
    return sum(@_)/@_;
}
