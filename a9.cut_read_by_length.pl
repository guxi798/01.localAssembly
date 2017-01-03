#!/usr/bin/perl -w
# run the script: perl 00.script/a9.cut_read_by_length.pl your_input_file your_output_file desired_length sequence_format

## read in parameters required by the script
my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;
my $read_len = shift @ARGV;
my $seq_format = lc(shift @ARGV);

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

if($seq_format eq "fastq"){
	while(my $line = <SRC>){
		chomp $line;
		if($line =~ /^@/){
			print TGT "$line\n";
			$line = <SRC>;
			chomp $line;
			my $subline = substr($line, 0, $read_len);
			print TGT "$subline\n";
			$line = <SRC>;
			print TGT "$line";
			$line = <SRC>;
			chomp $line;
			$subline = substr($line, 0, $read_len);
			print TGT "$subline\n";
		}
	}
}elsif($seq_format eq "fasta"){
	while(my $line = <SRC>){
		chomp $line;
		if($line =~ /^>/){
			print TGT "$line\n";
			$line = <SRC>;
			chomp $line;
			my $subline = substr($line, 0, $read_len);
			print TGT "$subline\n";
		}
	}
}else{
	die "Erorr: Please specify the read format: fastq or fasta";
}

close SRC;
close TGT;