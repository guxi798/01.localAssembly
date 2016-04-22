#!/usr/bin/perl -w
# run the script: time perl 00.script/04.retrieveunmap.reads.pl 3 1

use strict;

my $sam = shift @ARGV;
my $run = shift @ARGV;
my $long = shift @ARGV;	#1: for true, 0: for false
my $srcfolder = "04.retrieve.reads/03.bowtie.nucl/run.$run";
my $reffolder = "01.data/02.Fasta";

opendir(SRC, "$srcfolder");
my @subs = sort(grep(/^[0-9]+/, readdir(SRC)));
my @temp = @subs;
closedir SRC;

my @R = ("R1", "R2", "long");
my $i = 0;
my $prev = 0;

foreach my $sub (@subs){
	print ".Procesing $sub\n";
	my %hash =();
	$i ++;
	my $remove = 0;
	if($i > 2){$remove = shift @temp;}
	foreach my $R (@R){
		if(-e "$srcfolder/$sub/retrieved.$sub.$R.fasta" and -s "$srcfolder/$sub/retrieved.$sub.$R.fasta"){
			print "....Processing $R\n";
		## construct hash table for this group
			open(SUB, "$srcfolder/$sub/retrieved.$sub.$R.fasta");
			while(my $line = <SUB>){
				chomp $line;
				$line =~ s/\/1|\/2//g;
				if($line =~ /^>/ and not exists $hash{$line}){
					$hash{$line} = 0;
				}
			}
			close SUB;
			
		## go back to retrieved reads unmapped to this group
			print ".......Processing $sam\n";
			my $samfile = 0;
			my $tgtfile = 0;
			if($R eq "long" and $long == 1){
				if(!($sam =~ /F$|Fu$/)){	## need 454 data
					print "..........Skipping\n";
					next;
				}
				$samfile = "$reffolder/$sam/$sam.fna";
				$tgtfile = "$reffolder/$sam/$sam.tmp.$sub.fna";
			}elsif($R ne "long"){
				if($sam =~ /F$|Fu$/ and $long == 1){		## skip 454 data
					print "..........Skipping\n";
					next;
				}
				$tgtfile = "$reffolder/$sam/$sam.tmp.$sub.$R.fasta";
				if($i == 1){
					$samfile = "$reffolder/$sam/$sam.$R.fasta_pairs_$R.fasta";
				}else{
					$samfile = "$reffolder/$sam/$sam.tmp.$prev.$R.fasta";
					if($i > 2){ # remove tmp file two iterations earlier
						print ".............Removing $reffolder/$sam/$sam.tmp.$remove.$R.fasta\n";
						system("rm $reffolder/$sam/$sam.tmp.$remove.$R.fasta");
					}
				}
			}
			
			open(SAM, $samfile);
			open(TGT, ">$tgtfile");
			print ".............Writing $tgtfile\n";
			while(my $line = <SAM>){
				chomp $line;
				if($line =~ /^>/ and not exists $hash{$line}){
					print TGT "$line\n";
					$line = <SAM>;
					print TGT "$line";
				}
			}
			close SAM;
			close TGT;
			
			## retrieved unmapped single end reads
			if($R eq "R1"){
				if($i == 1){
					open(SAM, "$reffolder/$sam/$sam.$R.fasta_singles.fasta");
				}else{
					open(SAM, "$reffolder/$sam/$sam.tmp.$prev.single.fasta");
					if($i > 2){ # remove tmp file two iterations earlier
						print ".............Removing $reffolder/$sam/$sam.tmp.$remove.single.fasta\n";
						system("rm $reffolder/$sam/$sam.tmp.$remove.single.fasta");
					}
				}
				my $tgtfile = "$reffolder/$sam/$sam.tmp.$sub.single.fasta";
				open(TGT, ">$tgtfile");
				print ".............Writing $tgtfile\n";
				while(my $line = <SAM>){
					chomp $line;
					$line =~ s/\/1$|\/2$//;
					if($line =~ /^>/ and not exists $hash{$line}){
						print TGT "$line\n";
						$line = <SAM>;
						print TGT "$line";
					}
				}
				close SAM;
				close TGT;
			}
		}
	}
	$prev = $sub;
}

my $last = shift @temp;
my $cur = shift @temp;
print "last:$last\tcur:$cur\n";
system("cp $reffolder/$sam/$sam.tmp.$cur.* $srcfolder/99999/");
#system("rm $reffolder/$sam/$sam.tmp.$last.*");

#system("cat $srcfolder/99999/*.R1.* > $srcfolder/99999/retrieved.99999.R1.fasta");
#system("cat $srcfolder/99999/*.R2.* > $srcfolder/99999/retrieved.99999.R2.fasta");
#system("cat $srcfolder/99999/*.single.* >> $srcfolder/99999/retrieved.99999.R1.fasta");


