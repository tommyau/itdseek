#!/usr/bin/env perl
# Chun Hang AU (chau@hksh.com) (https://github.com/tommyau/itdseek)
# ITDsim - Simulate amplicon sequencing reads for FLT3 internal tandem duplication (FLT3 ITD)
# itdsim.pl - simulate wild-type and ITD sequences
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

############
# Amplicon and adapter definition
# The following is an example for the amplicon "FLT3.ITD.line.29.chr13.28607916.28608351_tile_2" in Illumina TruSight Myeloid Panel (reference: TruSight-Myeloid-Manifest.txt from Illumina)
# refseq is extracted using UCSC Genome Browser (>hg19_dna range=chr13:28607477-28609208 5'pad=0 3'pad=0 strand=+ repeatMasking=none)
my $setting_refseq = "GCAACTGTATCCCCAAATATCAAAGCCTATAAAAACTAAATAGGCCTGGCACAGTGGCTCACGCCTGTAATCCCAACACTTTGGGAGGCCGAGGCGGGCAGATAGTTTGAGCTCAGGAGTCCCAGATCAGCCTAGGCAACATGGTGAAACCCCGTCTCTACCAAAAATAAAAAACTTAGCTGAGCGTGGTGGTGCACGCCTGTAGCCCCAGCTGCTGAGGAGCCTGAGCCCAGGGGGTGGAGGCTGCAGTGAGCCATGATCACACTACTGTACTCCAGCCTAGGTGACAGAGTGAGACCCTGTCTCAAAAAAATAAAAGAAAATAAAAATAAACAAAGAGAGAAGTGGAAGAAGAGGTGGAGTTTTGTATTTATGACTTGAATTTTGTATTCATGACTGGGTTGACACCCCAATCCACTCCATTTTTAGCCTTGAAACATGGCAAACAGTAACCATTAAAAGGATGGAAAAGAGAAGAAGGCATGGGTGGGAAACTGTGCCTCCCATTTTTGTGCATCTTTGTTGCTGTCCTTCCACTATACTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCTTCAGAGATGAAATGATGAGTCAGTTAGGAATAGGCAGTTCTGCAGATAGAGGAAAGAATAATGAATTTTTACCTTTGCTTTTACCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACAAACACCAATTGTTGCATAGAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGCCTGCAACAAAAGAGTGTCACTCAGCGATGAAACAGAATTCCTGTGTGACATTATAAATAGTGGACAACTCATTATAATCTCTCACATCCTGTTTCAGTAATAATCATTTTCAGTCCTAACAACCACTCTACATATACTCTACTCCCCACAGACAATCAGGCAATGTCCCTGTAAAGGATACATTTCCTCCCTAGAAAATTGCGGATTATTCTCAATCCATTCTTTAAAACCATTTACTAGGGTAAATTTACAAGAATTACATCTGGTCCAGGCACGATGGCTCACGCCTGTAGTCCCAGCACTTTGGGAGGCCAAGATGGGAGGATCACTTGAGTCCAAGAATTAGACACCAGCCCAGGCAACACAGTGAAATCCCGTCTCTAAAAAAATTCAAAAATTAGCTGGGCGTGGTGGCAGGTGCCTGTAATCCCAGCTGCTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCTGGCAGAGGTGAGCCAAGATCACGCCACTGCACTGCAGCCTGGGTGACAGAGCGAGTCTCCATCTCAAAAAAAAAAATGCATGGATTTAAACCTACAGGATGTACAAGTCAGAGAAAAAAGAAGTTCATATCCCTTTTATGAAGGTCGCACATGAAAGTGGAGAAGTCAGTTTTGTGAGGTTGGTGATCCTAA";
my $setting_primerup = "TGCGTTCATCACTTTTCCAAAAGCACC";
my $primerup_start = 609;
my $primerup_end = 635;
my $setting_primerdown = "CACCTGTACCATCTGTAGCTGGCTTT";
my $primerdown_start = 837;
my $primerdown_end = 862;
my $setting_adapterup = "CAACGATCGTCGAAATTCGC";
my $setting_adapterdown = "AGATCGGAAGAGCGTCGTGTA";
############

# internal logic below
my $refseq = Bio::Seq->new(-display_id => "refseq", -seq => $setting_refseq);
my $SEQadapterup = Bio::Seq->new(-display_id => "adapterup", -seq => $setting_adapterup);
my $SEQprimerup = Bio::Seq->new(-display_id => "primerup", -seq => $setting_primerup);
my $SEQprimerdown = Bio::Seq->new(-display_id => "primerdown", -seq => $setting_primerdown);
my $SEQadapterdown = Bio::Seq->new(-display_id => "adapterdown", -seq => $setting_adapterdown);
# default parameters
my $R1readlen = -1;
my $R2readlen = -1;
my $ITDstartsstart = -201;
my $ITDstartsend = 201;
my $ITDlengthsstart = 1;
my $ITDlengthsend = 201;
my $debug = 0;

GetOptions ("r1=i" => \$R1readlen,
	    "r2=i" => \$R2readlen,
	    "itdstartsstart=i" => \$ITDstartsstart,
	    "itdstartsend=i" => \$ITDstartsend,
	    "itdlengthsstart=i" => \$ITDlengthsstart,
	    "itdlengthsend=i" => \$ITDlengthsend,
	    "debug" => \$debug,
	   );
my @ITDstarts = ($ITDstartsstart .. $ITDstartsend);
my @ITDlengths = ($ITDlengthsstart .. $ITDlengthsend);

ITDstart: foreach my $ITDstart (@ITDstarts) {
    my $internal_ITDstart;
    if ($ITDstart == 0){
	#die "possible ITDstart values: -inf .. -1, +1 .. inf; 0 position does not exist";
	next ITDstart; # to skip 0 to avoid panic
    } elsif ($ITDstart >= 1) {
	$internal_ITDstart = $ITDstart;
    } elsif ($ITDstart <= -1) {
	$internal_ITDstart = $ITDstart + 1;
    }
    ITDlength: foreach my $ITDlength (@ITDlengths) {
	my $ITDseq_start = $primerup_end + $internal_ITDstart;
	my $ITDseq_end = $ITDseq_start + $ITDlength - 1;
	
	my $preITDseq = $refseq->subseq(1, $ITDseq_start - 1);
	my $ITDunitseq = $refseq->subseq($ITDseq_start,$ITDseq_end);
	my $postITDseq = $refseq->subseq($ITDseq_end + 1, $refseq->length);
	
	# WT and ITD reference sequences
	my $WTseq =  $preITDseq . $ITDunitseq . $postITDseq;
	my $ITDseq =  $preITDseq . $ITDunitseq . $ITDunitseq . $postITDseq;
	# assert wild-type sequence is identical to reference sequence
	die "bug occured" unless $WTseq eq $refseq->seq;

	# simulate reads
	my ($WTseqR1, $WTseqR2) = &generate_amplicons($WTseq, $SEQprimerup->seq, $SEQprimerdown->seq);
	my ($ITDseqR1, $ITDseqR2) = &generate_amplicons($ITDseq, $SEQprimerup->seq, $SEQprimerdown->seq);
	die "unexpected case occured" if !defined $WTseqR1 || !defined $WTseqR2;
	next if !$debug && (!defined $ITDseqR1 || !defined $ITDseqR2);
	print join("\t", $ITDstart, $ITDlength, $preITDseq, $ITDunitseq, $postITDseq, $WTseq, $WTseqR1, $WTseqR2, $ITDseq, defined $ITDseqR1 ? $ITDseqR1 : "undef", defined $ITDseqR2 ? $ITDseqR2 : "undef") . "\n";
	
    }
}

sub generate_amplicons {
    my ($template, $upperprimer, $lowerprimer) = @_;
    my @PCRamplicons = $template =~ m/($upperprimer.+$lowerprimer)/g;
    if (scalar @PCRamplicons == 0) {
	return (undef, undef);
    } elsif (scalar @PCRamplicons == 1) {
	return &PCRamplicon_to_reads($PCRamplicons[0]);
    } else {
	die "unexpected case occured";
    }
}

sub PCRamplicon_to_reads {
    my ($PCRamplicon) = @_;
    # R2prelim = $seqprimerup . mutROI . $seqprimerdown. $seqadapterdown . "G"x1000
    my $seqR2prelim = Bio::Seq->new(-display_id => "R2prelim", -seq => $PCRamplicon . $SEQadapterdown->seq . "G" x 1000);
    # R2 = R2prelim subseq 1 ~~ R2readlen
    my $seqR2 = $R2readlen == -1 ? $seqR2prelim->seq : $seqR2prelim->subseq(1, $R2readlen);
    # R1prelim = "G"x1000 . $seqadapterup . $seqprimerup . mutROI . $seqprimerdown
    my $seqR1prelim = Bio::Seq->new(-display_id => "R1prelim", -seq =>  "G" x 1000 . $PCRamplicon);
    # R1 = R1 revcom subseq 1 ~~ R1readlen
    my $seqR1 = $R1readlen == -1 ? $seqR1prelim->revcom->seq : $seqR1prelim->revcom->subseq(1, $R1readlen);
    return ($seqR1, $seqR2);
}


exit;
