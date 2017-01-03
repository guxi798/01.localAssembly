#!/usr/bin/perl -w
# run the script: time perl AmericanDogwood/00.script/z3.readcount.len.depth.pl . readcount_summary.txt

use strict;

my $srcfolder = shift @ARGV;
my $tgtfile = shift @ARGV;

opendir(SRC, $srcfolder);
my @sp = sort(grep(/^[A-Z]/, readdir SRC));
closedir SRC;
open(TGT, ">$tgtfile");
print "$tgtfile\n";

foreach my $sp(@sp){
	# get sample name
	opendir(SRC, "$srcfolder/$sp/01.data/02.Fasta");
	my @subs = sort(grep(/\w+/, readdir SRC));
	my @temp = @subs;
	closedir SRC;
	
	# print header info
	print TGT "### $sp\n";
	print "### $sp\n";
	my @names1 = map {$_."_Mapped#"} @temp;
	my @names2 = map {$_."_Total#"} @temp;
	print TGT "contig_ID\tPtr_ID\tcontig_len\t",  join("\t", @names1), "\t", join("\t", @names2), "\n";
	#print "contig_ID\tPtr_ID\tcontig_len\t",  join("\t", @names1), "\t", join("\t", @names2), "\n";
	
	# get seq depth
	my @count1 = ();
	foreach my $sub(@subs){
		open(SUB, "$srcfolder/$sp/01.data/02.Fasta/$sub/$sub.R1.fasta");
		my $count = 0;
		while(my $line = <SUB>){
			if($line =~ /^>/){
				$count ++;
			}
		}
		
		push @count1, $count;
		close SUB;
	}
	
	# summarize read count
	@subs = @temp;
	my %hash  = ();
	foreach my $sub (@subs){
		open(SUB, "$srcfolder/$sp/13.abundance/01.SUT/$sub.count.txt");
		print "$srcfolder/$sp/13.abundance/01.SUT/$sub.count.txt\n";
		while(my $line = <SUB>){
			chomp $line;
			my @lines = split(/\s+/, $line);
			if(not exists $hash{$lines[0]}){
				$hash{$lines[0]} = [$lines[2]];
			}else{
				push @{$hash{$lines[0]}}, $lines[2];
			}
		}
	}
	
	# get contig length
	open(SRC, "$srcfolder/$sp/01.data/00.PriorData/contig2gene.txt");
	<SRC>;
	while(my $line = <SRC>){
		chomp $line;
		my @lines = split(/\s+/, $line);
		print TGT join("\t", @lines), "\t", join("\t", @{$hash{$lines[0]}}), "\t", join("\t", @count1), "\n";
	}
	close SRC;

}

close TGT;
