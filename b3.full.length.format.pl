#!/usr/bin/perl -w

use strict;

my $groupfile = shift @ARGV;
my $file = shift @ARGV;
my $out1 = shift @ARGV;
my $out2 = shift @ARGV;

open(SRC, $file);
open(GRP, $groupfile);
open(TGT, ">$out1");
open(TGT2, ">$out2");

my %group = ();
while(my $line = <GRP>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $group{$lines[0]}){
		$group{$lines[0]} = $lines[5];
	}
}

my %hash = ();
foreach my $line (<SRC>){
    chomp $line;
    my @lines = split(/\s+/, $line);
	my $qname = $lines[0];
	my $hitname = $lines[2];
	my $run = $lines[7] - 1;
	my $len = $lines[1];
	$len =~ s/len=//;
	
    if(not exists $hash{$hitname}){
        $hash{$hitname} = { run => $run, 
							query => [$qname],
							length => [$len]
						};
    }else{
		my $oldrun = $hash{$hitname}->{run};
        if($run == $oldrun){
            push @{$hash{$hitname}->{query}}, $qname;
			push @{$hash{$hitname}->{length}}, $len;
        }elsif($run > $oldrun){
			## obtain group info
			my $gr = $group{$hitname};
			
			## open old run blast file
			open(BLT1, "07.map.back/03.bowtie.nucl/run.$oldrun/$gr/$gr.contigs.blast.out");
			my $oldbits = 0;
			my $oldmismatch = 0;
			my $oldlen = 0;
			while(my $record = <BLT1>){
				if($record =~ /$hitname/){
					chomp $record;
					my @records = split(/\s+/, $record);
					my $bits = $records[11];
					my $mismatch = $records[4];
					if($bits > $oldbits){
						$oldbits = $bits;
						$oldmismatch = $mismatch;
						for(my $i=0; $i<scalar(@{$hash{$hitname}->{query}}); $i++){
							if(${$hash{$hitname}->{query}}[$i] eq $records[0]){
								$oldlen = ${$hash{$hitname}->{length}}[$i];
								last;
							}
						}
					}
				}
			}
			close BLT1;
			
			## open new run blast file
			open(BLT2, "07.map.back/03.bowtie.nucl/run.$run/$gr/$gr.contigs.blast.out");
			my $newbits = 0;
			my $newmismatch = 0;
			while(my $record = <BLT2>){
				if($record =~ /$hitname/){
					chomp $record;
					my @records = split(/\s+/, $record);
					my $bits = $records[11];
					my $mismatch = $records[4];
					if($bits > $newbits){
						$newbits = $bits;
						$newmismatch = $mismatch;
					}
				}
			}
			close BLT2;
			if($hitname eq "Potri.010G084700"){
				print "$hitname\t$oldrun\t$oldbits\t$oldlen\t$run\t$newbits\t$len\n";
			}
			if(($newbits > $oldbits) or ($newbits == $oldbits and $len >= $oldlen) or ($newbits >= $oldbits*0.95 and $newmismatch <= $oldmismatch and $len > $oldlen)){
				$hash{$hitname} = { run => $run, 
									query => [$qname],
									length => [$len]
								};
			}
        }
    }
}

foreach my $key (sort keys %hash){
    print TGT "$key\n";
    my $run = $hash{$key}->{run};
    foreach my $contig (@{$hash{$key}->{query}}){
        print TGT2 "$key\t$contig\t$run\n";
    }
}

close SRC;
close TGT;
close TGT2;

