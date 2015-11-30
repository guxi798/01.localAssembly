#!/usr/bin/perl -w
# run the script: 
# module load perl/5.20.2
# time perl 00.script/c7.get.blast.identity.pl 05.Full.Length/both.full.gene.txt 04.blastn/blast.xml.out ../03.Protein.Subset/07.map.back/02.blastn 05.Full.Length/both.full.identity.txt > 05.Full.Length/both.full.identity.out

use strict;
use Bio::SearchIO;
use List::Util qw(sum);

my $fullfile = shift @ARGV;
my $globlastfile = shift @ARGV;
my $locblastfolder = shift @ARGV;
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
#=pod
print "global blast.....\n";
my $globin = new Bio::SearchIO(-format 	=> 'blastxml',
								-file	=> $globlastfile);
my %globblast = ();
while(my $result = $globin->next_result){
	my @names = split(/\s+/, $result->query_description);
	my $len = $names[1];
	$len =~ s/len=//;
	if(not exists $globfull{$names[0]}){next;}
	while(my $hit = $result->next_hit){
		my @hitname = split(/\./, $hit->name);
		#pop @hitname;
		my $hitname = join(".", @hitname);
		print "$names[0]\t$len\t$hitname\t";
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		print $identity, "\t",$identity/$len,"\n";
		if(not exists $globblast{$hitname}){
			$globblast{$hitname}{'I'} = [$identity/$len];
			$globblast{$hitname}{'L'} = [$len];
		}else{
			push @{$globblast{$hitname}{'I'}}, $identity/$len;
			push @{$globblast{$hitname}{'L'}}, $len;
		}
	}
}
#=cut

## get in local blast record
print "local blast.....\n";
my %locblast = ();
opendir(SRC, $locblastfolder);
my @runs = sort grep(/run/, readdir SRC);
foreach my $i (0..(@runs-1)){
	$runs[$i] =~ s/run\.//;
}
closedir SRC;
@runs = sort {$a <=> $b} @runs;

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
			if(not exists $locfull{$fullrun}{$group}{$names[0]}){next;}
			my $len = $names[1];
			$len =~ s/len=//;
			while(my $hit = $result->next_hit){
				my $hitname = $hit->name;
				print "$names[0]\t$len\t$hitname\t";
				my $identity = scalar($hit->seq_inds('query', 'identical'));
				print $identity, "\t", $identity/$len,"\n";
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

#=pod
## re-read fully recovered genes 
open(TGT, ">$tgtfile");
foreach my $key (sort(keys %genefull)){
	my $globident = 0;
	my $globlen = 0;
	my $locident = 0;
	my $loclen = 0;
	
	if(exists $globblast{$key}){
		my @glob = @{$globblast{$key}{'I'}};
		my @len = @{$globblast{$key}{'L'}};
		$globident = &mean(@glob);
		$globlen = &mean(@len);
	}
	
    if(exists $locblast{$key}){
        my @temp = ();
        my @len = ();
        foreach my $contig (sort(keys %{$locblast{$key}})){
            my @loc = @{$locblast{$key}{$contig}{'I'}};
            push @temp, &mean(@loc);
            push @len, $locblast{$key}{$contig}{'L'};
        }
        $locident = &mean(@temp);
        $loclen = &mean(@len);
    }
	
	print TGT "$key\t$globlen\t$globident\t$loclen\t$locident\n";
}
#=cut

############################
sub mean {
    return sum(@_)/@_;
}
