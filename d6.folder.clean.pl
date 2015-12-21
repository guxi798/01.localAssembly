#!/usr/bin/perl -w
use strict;

my $srcfolder = shift @ARGV;

opendir(SRC, $srcfolder);
my @subs = sort(grep(/^[A-Z]\w.+/, readdir(SRC)));
my @R = ("R1", "R2");

system("rm -rf 00.script/d6.shell.script");
system("mkdir -p 00.script/d6.shell.script");

foreach my $sub(@subs)
{
    foreach my $R (@R){
        if($sub =~ /F$|Fu$|R$|OrAe1G$/){next;}
        #if($sub =~ /F$|Fu$|R$/){next;}
        my $shell = "00.script/d6.shell.script/d6.final.clean.$sub.$R.sh";

        open(SHL,">$shell") or die "Cannot write $shell: $!";
        print SHL "#PBS -S /bin/bash\n";
        print SHL "#PBS -q batch\n";
        print SHL "#PBS -N d6.final.clean.$sub.$R.sh\n";
        print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
        print SHL "#PBS -l walltime=48:00:00\n";
        print SHL "#PBS -l mem=60gb\n";
        print SHL "\n";
        print SHL "cd \$PBS_O_WORKDIR\n";

        print SHL "time perl 00.script/d6.final.clean.pl $reffolder/$sub/$sub.$R.fastq $srcfolder/$sub/$sub.$R.fastq $srcfolder/$sub/$sub.$R.clean\n";
        close(SHL);

        system("chmod u+x $shell");
        system("qsub $shell");
}
}

close(SRC);
