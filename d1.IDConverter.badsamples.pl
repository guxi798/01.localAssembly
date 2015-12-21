#!/usr/bin/perl -w

use strict;

my $srcfile = shift @ARGV;
my $tgtfile = shift @ARGV;

open(SRC, $srcfile);
open(TGT, ">$tgtfile");

while(my $line = <SRC>)

close SRC;
close TGT;