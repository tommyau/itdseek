#!/usr/bin/env perl
# Chun Hang AU (chau@hksh.com) (https://github.com/tommyau/itdseek)
# ITDsim - Simulate amplicon sequencing reads for FLT3 internal tandem duplication (FLT3 ITD)
# itdsim2fastq.pl - convert itdsim.pl output to FASTQ
use strict;
use warnings;
use IO::File;
use Getopt::Long;
my $r1 = 275;
my $r2 = 275;
my $depthwt = 1000;
my $depthitd = 1000;
GetOptions ("r1=i" => \$r1,
	    "r2=i" => \$r2,
	    "depthwt=i" => \$depthwt, 
	    "depthitd=i" => \$depthitd, 
	   );

while (<>) {
    chomp;
    my @f=split("\t");
    my $fhR1 = IO::File->new(">ITD.$f[0].$f[1].R1.fastq"); 
    my $fhR2 = IO::File->new(">ITD.$f[0].$f[1].R2.fastq"); 
    my $wtr1seq = substr($f[6],0,$r1);
    my $wtr2seq = substr($f[7],0,$r2);
    my $itdr1seq = substr($f[9],0,$r1);
    my $itdr2seq = substr($f[10],0,$r2);
    foreach my $i (1..$depthwt) {
        my $identifier = "$f[0]:$f[1]:$i";
	print {$fhR1} "\@WT:$identifier\n$wtr1seq\n+\n".("I" x length($wtr1seq))."\n";
	print {$fhR2} "\@WT:$identifier\n$wtr2seq\n+\n".("I" x length($wtr2seq))."\n";
    }
    foreach my $i (1..$depthitd) {
        my $identifier = "$f[0]:$f[1]:$i";
	print {$fhR1} "\@ITD:$identifier\n$itdr1seq\n+\n".("I" x length($itdr1seq))."\n";
	print {$fhR2} "\@ITD:$identifier\n$itdr2seq\n+\n".("I" x length($itdr2seq))."\n";
    }
}
