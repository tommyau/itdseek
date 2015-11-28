#!/usr/bin/env perl
# Chun Hang AU (chau@hksh.com) (https://github.com/tommyau/itdseek)
# ITDseek - Detect FLT3 internal tandem duplication (FLT3 ITD) in amplicon sequencing reads
# Example usage:
# /home/adminrig/tools/samtools-0.1.19/samtools view -f 0x0100 1472.bwa-mem.aln-pe.L5.bam chr13:28607161-28609590 | perl itdseek.pl --refseq /home/adminrig/tools/GenomeAnalysisTK-2.8-1-g932cd3a/bundle_2.8_hg19/ucsc.hg19.fasta --samtools /home/adminrig/tools/samtools-0.1.19/samtools > 1472.itdseek.vcf
use strict;
use warnings;
use Getopt::Long;

# configuration for default paths to reference sequence and samtools
my $default_refseq = "/home/adminrig/tools/GenomeAnalysisTK-2.8-1-g932cd3a/bundle_2.8_hg19/ucsc.hg19.fasta";
my $default_samtools = "/home/adminrig/tools/samtools-0.1.19/samtools";

##
# internal logic below
my %variants;
my ($refseq, $samtools);
GetOptions ("refseq=s" => \$refseq,
	    "samtools=s" => \$samtools,);
$refseq = $default_refseq if !defined $samtools;
$samtools = $default_samtools if !defined $samtools;

ALIGNMENT:while(<>) {
    chomp;
    my @fields = split "\t", $_;
    # skip SAM header lines starting with @, in case user pipes in whole raw SAM file
    next ALIGNMENT if substr($fields[0],0,1) eq "@";
    # NOT_PRIMARY mask is 0x0100 (http://cpansearch.perl.org/src/LDS/Bio-SamTools-1.41/lib/Bio/DB/Sam/Constants.pm)
    next ALIGNMENT unless $fields[1] & 0x0100;
    my $cigar = $fields[5];
    my $readname = $fields[0];
    my $readseq = $fields[9];
    my $refseq = $fields[2];
    #REVERSED mask is 0x0010
    my $direction = $fields[1] & 0x0010 ? "-" : "+";
    my $start = $fields[3];

    my (@sa, $sa_tag);
    OPTIONALFIELD:foreach my $i (11..$#fields) {
	if ($fields[$i]=~/^SA:Z:(.*)$/) {
	    @sa = split /\;/, $1;
	    $sa_tag=$fields[$i];
	    last OPTIONALFIELD;
	}
    }
    # expect only one secondary alignment in SA
    next ALIGNMENT unless scalar @sa == 1;
    my @sa0_fields = split /\,/, $sa[0];
    # expect same chromosome as primary alignment
    next ALIGNMENT unless $refseq eq $sa0_fields[0];
    # expect same strand as primary alignment
    next ALIGNMENT unless $direction eq $sa0_fields[2];
    
    
    my ($clip_5prime, $clip_3prime, $clip_type) = &get_clipping_details($start, $cigar);
    my ($pm_clip_5prime, $pm_clip_3prime, $pm_clip_type) = &get_clipping_details($sa0_fields[1], $sa0_fields[3]);
    next ALIGNMENT if !defined $clip_5prime || !defined $pm_clip_5prime;
    
    my $itd_length = abs($pm_clip_3prime - $clip_3prime);
    my ($seqstart, $seqstop) = ($pm_clip_5prime < $clip_3prime) ? ($pm_clip_5prime, $clip_3prime) : ($clip_3prime, $pm_clip_5prime);
    my $sequencestring = sprintf("%s:%d-%d", $refseq, $seqstart, $seqstop);
    
    #print join("\t", "#".$refseq, $clip_5prime, $clip_3prime, $direction, $clip_type, $readname, $cigar, $sa_tag, $pm_clip_5prime, $pm_clip_3prime, $pm_clip_type, $itd_length, $seqstart, $seqstop)."\n";

    my $variantid=sprintf("%s:%d-%d", $refseq,$seqstart,$seqstop);
    if (!defined $variants{$variantid}){
       $variants{$variantid} = [0, 0, $refseq,$seqstart,$seqstop];
    }
    if ($direction eq "+"){
       $variants{$variantid}->[0]++;
    } elsif ($direction eq "-"){
       $variants{$variantid}->[1]++;
    }

}

# VCF output
print "##fileformat=VCFv4.1\n";
print "##source=ITDseek\n";
print "##reference=file://$refseq\n";
print "##INFO=<ID=DP2,Number=2,Type=Integer,Description=\"# alt-foward and alt-reverse reads\">\n";
print "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"length of ITD\">\n";
print "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"sequence of ITD\">\n";
print join("\t", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")."\n";
foreach my $variantid (sort {$variants{$b}->[0]+$variants{$b}->[1] <=> $variants{$a}->[0]+$variants{$a}->[1]} keys %variants) {
    my $seq = &faidx($variantid);
    my $ref = substr($seq,-1);
    my $alt = $ref.$seq;
    print join("\t", $variants{$variantid}->[2], $variants{$variantid}->[4], ".", $ref, $alt, $variants{$variantid}->[0] + $variants{$variantid}->[1], ".", sprintf("DP2=%d,%d;LEN=%d;SEQ=%s", $variants{$variantid}->[0], $variants{$variantid}->[1], $variants{$variantid}->[4] - $variants{$variantid}->[3], $seq))."\n";
}



sub get_clipping_details {
    my ($alignment_start, $cigar) = @_;
    my ($clip_5prime, $clip_3prime, $clip_type);
    if ($cigar =~ m/^\d+([SH])\d+M/){
	$clip_5prime = $alignment_start - 1;
	$clip_3prime = $clip_5prime + 1;
	$clip_type = $1;
    } elsif ($cigar =~ m/(\d+)M\d+([SH])$/){
	$clip_5prime = $alignment_start + $1 - 1;
	$clip_3prime = $clip_5prime + 1;
	$clip_type = $2;
    }
    return ($clip_5prime, $clip_3prime, $clip_type);
}

sub faidx {
    my ($region) = @_;
    my $output = `$samtools faidx $refseq $region 2>&1`;
    my @outputlines = split("\n",$output);
    die "Error: unexpected output from samtools faidx" unless $outputlines[0] eq ">$region";
    return join("",@outputlines[1..$#outputlines]);
}