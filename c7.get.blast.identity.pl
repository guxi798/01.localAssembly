#!/usr/bin/perl -w
# run the script: 
# module load perl/5.20.2
# time perl 00.script/c7.get.blast.identity.pl 08.GlobLocal.Comparison 05.Trinity/blastx.xml.out 07.map.back/03.bowtie.nucl > 08.GlobLocal.Comparison/both.full.identity.out

use strict;
use Bio::SearchIO;
use List::Util qw(sum);

my $outfolder = shift @ARGV;
my $globlastfile = shift @ARGV;
my $locblastfolder = shift @ARGV;
my $bothfull = "$outfolder/both.full.gene.txt";
my $globfull = "$outfolder/global.uniq.gene.txt";
my $locfull = "$outfolder/local.uniq.gene.txt";
my $tgtfile = "$outfolder/both.full.identity.txt";
my $globtgt = "$outfolder/global.uniq.identity.txt";
my $loctgt = "$outfolder/local.uniq.identity.txt";

## put fully recovered genes into hash table
open(FUL, $bothfull);<FUL>;
open(FUL1, $globfull);<FUL1>;
open(FUL2, $locfull);<FUL2>;

my %globfull = ();
my %locfull = ();
my %genefull = ();
my %globuniq = ();
my %locuniq = ();

foreach my $line(<FUL>){
	chomp $line;
	my @lines = split(/\t/, $line);
	
	## store unique gene ID into hash table
	if(not exists $genefull{$lines[0]}){
		$genefull{$lines[0]} = 0;
	}else{
		print "This gene appeared more than once: NEED CHECK $lines[0]!\n";
	}
	
	## store global contig ID
	my @glob = split(/,/, $lines[1]);
	foreach my $glob (@glob){
		if(not exists $globfull{$glob}){
			$globfull{$glob} = 0;
		}
	}
	
	## store local contig ID, as diff runs/groups can share the same ID, need to 
	## include both run and group info for unique contig ID
	my @loc = split(/,/, $lines[2]);
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

foreach my $line(<FUL1>){
	chomp $line;
	my @lines = split(/\t/, $line);
	## store unique gene ID into hash table
	if(not exists $genefull{$lines[0]}){
		$genefull{$lines[0]} = 0;
	}else{
		print "This gene appeared more than once: NEED CHECK $lines[0]!\n";
	}

	## store global contig ID
	my @glob = split(/,/, $lines[1]);
	foreach my $glob (@glob){
		if(not exists $globuniq{$glob}){
			$globuniq{$glob} = 0;
		}
	}
}
close FUL1;

foreach my $line(<FUL2>){
	chomp $line;
	my @lines = split(/\t/, $line);
	## store unique gene ID into hash table
	if(not exists $genefull{$lines[0]}){
		$genefull{$lines[0]} = 0;
	}else{
		print "This gene appeared more than once: NEED CHECK $lines[0]!\n";
	}

	## store local contig ID, as diff runs/groups can share the same ID, need to 
	## include both run and group info for unique contig ID
	my @loc = split(/,/, $lines[1]);
	foreach my $temp (@loc){
		my @temp = split(/:/, $temp);
		my $loc = shift @temp;
		my $group = shift @temp;
		my $run = shift @temp;
		
		if(not exists $locuniq{$run}{$group}{$loc}){
			$locuniq{$run}{$group}{$loc} = 0;
		}
	}
}
close FUL2;

## get in global blast record
#=pod
print "global blast.....\n";
my $globin = new Bio::SearchIO(-format 	=> 'blastxml',
								-file	=> $globlastfile);
my %globblast = ();
my %globblastuniq = ();
while(my $result = $globin->next_result){
	my @names = split(/\s+/, $result->query_description);
	my $len = $result->query_length;
	if(not exists $globfull{$names[0]} and not exists $globuniq{$names[0]}){next;}
	while(my $hit = $result->next_hit){
		my $hitname = $hit->name;
		print "$names[0]\t$len\t$hitname\t";
		my $identity = scalar($hit->seq_inds('query', 'identical'));
		my $conserved = scalar($hit->seq_inds('query', 'conserved'));
		my $frac_iden = $hit->frac_identical;
		my $frac_cons = $hit->frac_conserved;
		print "$identity\t", $identity/$len, "\t$frac_iden\t$conserved\t", $conserved/$len, "\t$frac_cons\n";
		if(exists $globfull{$names[0]}){
			if(not exists $globblast{$hitname}){
				$globblast{$hitname}{'I'} = [$identity/$len];
				$globblast{$hitname}{'IE'} = [$frac_iden];
				$globblast{$hitname}{'C'} = [$conserved/$len];
				$globblast{$hitname}{'CE'} = [$frac_cons];
				$globblast{$hitname}{'L'} = [$len];
			}else{
				push @{$globblast{$hitname}{'I'}}, $identity/$len;
				push @{$globblast{$hitname}{'IE'}}, $frac_iden;
				push @{$globblast{$hitname}{'C'}}, $conserved/$len;
				push @{$globblast{$hitname}{'CE'}}, $frac_cons;
				push @{$globblast{$hitname}{'L'}}, $len;
			}
		}elsif(exists $globuniq{$names[0]}){
			if(not exists $globblastuniq{$hitname}){
				$globblastuniq{$hitname}{'I'} = [$identity/$len];
				$globblastuniq{$hitname}{'IE'} = [$frac_iden];
				$globblastuniq{$hitname}{'C'} = [$conserved/$len];
				$globblastuniq{$hitname}{'CE'} = [$frac_cons];
				$globblastuniq{$hitname}{'L'} = [$len];
			}else{
				push @{$globblastuniq{$hitname}{'I'}}, $identity/$len;
				push @{$globblastuniq{$hitname}{'IE'}}, $frac_iden;
				push @{$globblastuniq{$hitname}{'C'}}, $conserved/$len;
				push @{$globblastuniq{$hitname}{'CE'}}, $frac_cons;
				push @{$globblastuniq{$hitname}{'L'}}, $len;
			}
		}
	}
}
#=cut

## get in local blast record
print "local blast.....\n";
my %locblast = ();
my %locblastuniq = ();

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
		my $fullrun = $run + 1;
		if(not exists $locfull{$fullrun}{$group} and not exists $locuniq{$fullrun}{$group}){next;}	# save some computational time
		print "run.$run\t$group\n";
		my $locin = new Bio::SearchIO(-format => 'blastxml', 
				-file => "$locblastfolder/run.$run/$group/$group.contigs.blast.xml.out");
		while(my $result = $locin->next_result){
			my @names = split(/\s+|\|/, $result->query_description);
			if(not exists $locfull{$fullrun}{$group}{$names[0]} and not exists $locuniq{$fullrun}{$group}{$names[0]}){next;}
			my $len = $result->query_length;
			while(my $hit = $result->next_hit){
				my $hitname = $hit->name;
				print "$names[0]\t$len\t$hitname\t";
				my $identity = scalar($hit->seq_inds('query', 'identical'));
				my $conserved = scalar($hit->seq_inds('query', 'conserved'));
				my $frac_iden = $hit->frac_identical;
				my $frac_cons = $hit->frac_conserved;
				print "$identity\t", $identity/$len, "\t$frac_iden\t$conserved\t", $conserved/$len, "\t$frac_cons\n";
				if(exists $locfull{$fullrun}{$group} and exists $locfull{$fullrun}{$group}{$names[0]}){
					if(not exists $locblast{$hitname}{$names[0]}){
						$locblast{$hitname}{$names[0]}{'I'} = [$identity/$len];
						$locblast{$hitname}{$names[0]}{'IE'} = [$frac_iden];
						$locblast{$hitname}{$names[0]}{'C'} = [$conserved/$len];
						$locblast{$hitname}{$names[0]}{'CE'} = [$frac_cons];
						$locblast{$hitname}{$names[0]}{'L'} = [$len];
					}else{
						push @{$locblast{$hitname}{$names[0]}{'I'}}, $identity/$len;
						push @{$locblast{$hitname}{$names[0]}{'IE'}}, $frac_iden;
						push @{$locblast{$hitname}{$names[0]}{'C'}}, $conserved/$len;
						push @{$locblast{$hitname}{$names[0]}{'CE'}}, $frac_cons;
						push @{$locblast{$hitname}{$names[0]}{'L'}}, $len;
						print "Warning: $hitname + $names[0] has occurred before\n";
					}
				}elsif(exists $locuniq{$fullrun}{$group} and exists $locuniq{$fullrun}{$group}{$names[0]}){
					if(not exists $locblastuniq{$hitname}{$names[0]}){
						$locblastuniq{$hitname}{$names[0]}{'I'} = [$identity/$len];
						$locblastuniq{$hitname}{$names[0]}{'IE'} = [$frac_iden];
						$locblastuniq{$hitname}{$names[0]}{'C'} = [$conserved/$len];
						$locblastuniq{$hitname}{$names[0]}{'CE'} = [$frac_cons];
						$locblastuniq{$hitname}{$names[0]}{'L'} = [$len];
					}else{
						push @{$locblastuniq{$hitname}{$names[0]}{'I'}}, $identity/$len;
						push @{$locblastuniq{$hitname}{$names[0]}{'IE'}}, $frac_iden;
						push @{$locblastuniq{$hitname}{$names[0]}{'C'}}, $conserved/$len;
						push @{$locblastuniq{$hitname}{$names[0]}{'CE'}}, $frac_cons;
						push @{$locblastuniq{$hitname}{$names[0]}{'L'}}, $len;
						print "Warning: $hitname + $names[0] has occurred before\n";
					}
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
print TGT "Gene\tGlob_len\tGlob_ident\tGlob_identE\tGlob_cons\tGlob_consE\tLoc_len\tLoc_ident\tLoc_identE\tLoc_cons\tLoc_consE\n";
open(TGT1, ">$globtgt");
print TGT1 "Gene\tGlob_len\tGlob_ident\tGlob_identE\tGlob_cons\tGlob_consE\n";
open(TGT2, ">$loctgt");
print TGT2 "Gene\tLoc_len\tLoc_ident\tLoc_identE\tLoc_cons\tLoc_consE\n";

foreach my $key (sort(keys %genefull)){
	my ($globident, $globidentE, $globcons, $globconsE, $globlen) = (0);
	my ($locident,  $locidentE, $loccons, $locconsE, $loclen) = (0);

	if(exists $globblast{$key} and exists $locblast{$key}){
		$globident = &mean(@{$globblast{$key}{'I'}});
		$globidentE = &mean(@{$globblast{$key}{'IE'}});
		$globcons = &mean(@{$globblast{$key}{'C'}});
		$globconsE = &mean(@{$globblast{$key}{'CE'}});
		$globlen = &mean(@{$globblast{$key}{'L'}});

        my (@tempI, @tempIE, @tempC, @tempCE, @len) = (());
        foreach my $contig (sort(keys %{$locblast{$key}})){
            push @tempI, &mean(@{$locblast{$key}{$contig}{'I'}});
			push @tempIE, &mean(@{$locblast{$key}{$contig}{'IE'}});
			push @tempC, &mean(@{$locblast{$key}{$contig}{'C'}});
			push @tempCE, &mean(@{$locblast{$key}{$contig}{'CE'}});
            push @len, &mean(@{$locblast{$key}{$contig}{'L'}});
        }
        $locident = &mean(@tempI);
		$locidentE = &mean(@tempIE);
		$loccons = &mean(@tempC);
		$locconsE = &mean(@tempCE);
        $loclen = &mean(@len);
		
		print TGT "$key\t$globlen\t$globident\t$globidentE\t$globcons\t$globconsE\t$loclen\t$locident\t$locidentE\t$loccons\t$locconsE\n";
    }

	## get info for unique genes
	if(exists $globblastuniq{$key}){
		$globident = &mean(@{$globblastuniq{$key}{'I'}});
		$globidentE = &mean(@{$globblastuniq{$key}{'IE'}});
		$globcons = &mean(@{$globblastuniq{$key}{'C'}});
		$globconsE = &mean(@{$globblastuniq{$key}{'CE'}});
		$globlen = &mean(@{$globblastuniq{$key}{'L'}});
		
		print TGT1 "$key\t$globlen\t$globident\t$globidentE\t$globcons\t$globconsE\n";
	}

    if(exists $locblastuniq{$key}){
        my (@tempI, @tempIE, @tempC, @tempCE, @len) = (());
        foreach my $contig (sort(keys %{$locblastuniq{$key}})){
            push @tempI, &mean(@{$locblastuniq{$key}{$contig}{'I'}});
			push @tempIE, &mean(@{$locblastuniq{$key}{$contig}{'IE'}});
			push @tempC, &mean(@{$locblastuniq{$key}{$contig}{'C'}});
			push @tempCE, &mean(@{$locblastuniq{$key}{$contig}{'CE'}});
            push @len, &mean(@{$locblastuniq{$key}{$contig}{'L'}});
        }
        $locident = &mean(@tempI);
		$locidentE = &mean(@tempIE);
		$loccons = &mean(@tempC);
		$locconsE = &mean(@tempCE);
        $loclen = &mean(@len);
		
		print TGT2 "$key\t$loclen\t$locident\t$locidentE\t$loccons\t$locconsE\n";
    }
}

close TGT;
close TGT1;
close TGT2;

############################
sub mean {
    my $mean = sum(@_)/@_;
	$mean = sprintf("%.4f", $mean);
}
