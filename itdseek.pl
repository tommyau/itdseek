#!/usr/bin/env perl
# Chun Hang AU (chau@hksh.com) (https://github.com/tommyau/itdseek)
# ITDseek - Detect FLT3 internal tandem duplication (FLT3 ITD) in amplicon sequencing reads
# Example usage:
# /home/adminrig/tools/samtools-1.3.1/samtools view 1472.bwa-mem.aln-pe.L5.bam chr13:28607161-28609590 | perl itdseek.pl --refseq /home/adminrig/tools/GenomeAnalysisTK-2.8-1-g932cd3a/bundle_2.8_hg19/ucsc.hg19.fasta --samtools /home/adminrig/tools/samtools-1.3.1/samtools --bam 1472.bwa-mem.aln-pe.L5.bam > 1472.itdseek.vcf
use strict;
use warnings;
use Getopt::Long;

# configuration for default paths to reference sequence and samtools
my $default_refseq = "/home/adminrig/tools/GenomeAnalysisTK-2.8-1-g932cd3a/bundle_2.8_hg19/ucsc.hg19.fasta";
my $default_samtools = "/home/adminrig/tools/samtools-1.3.1/samtools";
my $max_samtools_depth = 500000; # overriding samtools depth default cap of 8000x depth

##
# internal logic below
my $bam;
my %variants; # from clipping mode
my %insertionvariants; # from insertion mode
my $min_insertion_len = 3;
my ($refseq, $samtools);
GetOptions ("refseq=s" => \$refseq,
	    "samtools=s" => \$samtools,
	    "bam=s" => \$bam,);
$refseq = $default_refseq if !defined $samtools;
$samtools = $default_samtools if !defined $samtools;

ALIGNMENT:while(<>) {
    chomp;
    my @fields = split "\t", $_;
    # skip SAM header lines starting with @, in case user pipes in whole raw SAM file
    next ALIGNMENT if substr($fields[0],0,1) eq "@";
    # NOT_PRIMARY mask is 0x0100 (http://cpansearch.perl.org/src/LDS/Bio-SamTools-1.41/lib/Bio/DB/Sam/Constants.pm)
    if ($fields[1] & 0x0100) {
	# clipping mode - supplementary alignments
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
	
	my $variantid=sprintf("%s:%d-%d", $refseq,$seqstart,$seqstop);
	if (!defined $variants{$variantid}){
	   $variants{$variantid} = [0, 0, $refseq,$seqstart,$seqstop];
	}
	if ($direction eq "+"){
	   $variants{$variantid}->[0]++;
	} elsif ($direction eq "-"){
	   $variants{$variantid}->[1]++;
	}
    } else {
	# insertion mode - plain CIGAR insertion in main alignments
	if ($fields[5] =~ /I/) {
	    # insertion found in CIGAR
	    my $reconstructed_alignment = &reconstruct_alignment(\@fields);
	    # loop through each CIGAR operations
	    my $pointer = 1;
	    CIGAROP:for (my $i = 0; $i<=$#{$reconstructed_alignment->{cigar}}; $i++) {
		my $start = $pointer;
		$pointer += $reconstructed_alignment->{cigar}->[$i]->[0];
		next CIGAROP unless $reconstructed_alignment->{cigar}->[$i]->[1] eq "I" && $reconstructed_alignment->{cigar}->[$i]->[0] >= $min_insertion_len;
		my $end = $pointer - 1;
		my $insertion_seq = join("", @{$reconstructed_alignment->{cigar_readseq}}[$start..$end]);
		# the insertion is between $reconstructed_alignment->{cigar_refpos}->[$start-1] and $reconstructed_alignment->{cigar_refpos}->[$end+1]
		if ($reconstructed_alignment->{cigar_refpos}->[$start-1] eq "*") {
		    warn "WARNING: unknown reference position (it is soft-clipped?): ".sprintf("\$reconstructed_alignment->{cigar_refpos}->[\$start-1]=%s; \$reconstructed_alignment->{cigar}->[\$start-1]=%s", $reconstructed_alignment->{cigar_refpos}->[$start-1], $reconstructed_alignment->{cigar}->[$start-1]);
		    next CIGAROP;
		}
		my $variantid = sprintf("%s:%d:%s", $fields[2], $reconstructed_alignment->{cigar_refpos}->[$start-1], $insertion_seq);
		if (!defined $insertionvariants{$variantid}){
		   $insertionvariants{$variantid} = [0, 0, $fields[2], $reconstructed_alignment->{cigar_refpos}->[$start-1], $insertion_seq];
		}
		if ($reconstructed_alignment->{direction} eq "+"){
		   $insertionvariants{$variantid}->[0]++;
		} elsif ($reconstructed_alignment->{direction} eq "-"){
		   $insertionvariants{$variantid}->[1]++;
		}
	    }
	    
	} else {
	    next ALIGNMENT;
	}
    }
}
# VCF output
print "##fileformat=VCFv4.1\n";
print "##source=ITDseekV1.2\n";
print "##reference=file://$refseq\n";
print "##INFO=<ID=DP2,Number=2,Type=Integer,Description=\"# alt-foward and alt-reverse reads\">\n";
print "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"length of ITD\">\n";
print "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"sequence of ITD\">\n";
print "##INFO=<ID=CLIPPING,Number=0,Type=Flag,Description=\"ITD is detected as soft-clipping\">\n";
print "##INFO=<ID=INSERTION,Number=0,Type=Flag,Description=\"ITD is detected as insertion\">\n";
print "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"ITD allele fraction\">\n";
print join("\t", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")."\n";
# clipping mode
VAR:foreach my $variantid (sort {$variants{$b}->[0]+$variants{$b}->[1] <=> $variants{$a}->[0]+$variants{$a}->[1]} keys %variants) {
    my $seq = &faidx($variantid);
    my $ref = substr($seq,-1);
    my $alt = $ref.$seq;
    my $region = sprintf("%s:%d-%d", $variants{$variantid}->[2], $variants{$variantid}->[4] - 0, $variants{$variantid}->[4] + 0);
    my $region_depth = depth($region);
    if ($region_depth <= 0) {
	warn "WARNING: unexpected region depth $region_depth for $region, the variant is skipped (id: $variantid)";
	next VAR;
    }
    my $vaf = ($variants{$variantid}->[0] + $variants{$variantid}->[1]) / $region_depth;
    print join("\t", $variants{$variantid}->[2], $variants{$variantid}->[4], ".", $ref, $alt, $variants{$variantid}->[0] + $variants{$variantid}->[1], ".", sprintf("DP2=%d,%d;LEN=%d;SEQ=%s;CLIPPING;VAF=%0.2f", $variants{$variantid}->[0], $variants{$variantid}->[1], length($seq), $seq, $vaf))."\n";
}
# insertion mode
VAR:foreach my $variantid (sort {$insertionvariants{$b}->[0]+$insertionvariants{$b}->[1] <=> $insertionvariants{$a}->[0]+$insertionvariants{$a}->[1]} keys %insertionvariants) {
    my $ref = &faidx(sprintf("%s:%d-%d", @{$insertionvariants{$variantid}}[2,3,3]));
    my $alt = $ref.$insertionvariants{$variantid}->[4];
    my $region = sprintf("%s:%d-%d", $insertionvariants{$variantid}->[2], $insertionvariants{$variantid}->[3] - 0, $insertionvariants{$variantid}->[3] + 0);
    my $region_depth = depth($region);
    if ($region_depth <= 0) {
	warn "WARNING: unexpected region depth $region_depth for $region, the variant is skipped (id: $variantid)";
	next VAR;
    }
    my $vaf = ($insertionvariants{$variantid}->[0] + $insertionvariants{$variantid}->[1]) / $region_depth;
    print join("\t", $insertionvariants{$variantid}->[2], $insertionvariants{$variantid}->[3], ".", $ref, $alt, $insertionvariants{$variantid}->[0] + $insertionvariants{$variantid}->[1], ".", sprintf("DP2=%d,%d;LEN=%d;SEQ=%s;INSERTION;VAF=%0.2f", $insertionvariants{$variantid}->[0], $insertionvariants{$variantid}->[1], length($insertionvariants{$variantid}->[4]), $insertionvariants{$variantid}->[4], $vaf))."\n";
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

my %faidxcache;
my $faidxcache_count = 0;
sub faidx {
    my ($region) = @_;
    return $faidxcache{$region} if exists $faidxcache{$region};
    my $output = `$samtools faidx $refseq $region 2>&1`;
    my @outputlines = split("\n",$output);
    die "Error: unexpected output from samtools faidx" unless $outputlines[0] eq ">$region";
    my $seq = join("",@outputlines[1..$#outputlines]);
    $faidxcache{$region} = $seq;
    $faidxcache_count++;
    return $seq;
}

sub parse_cigar {
    my ($cigar_string) = @_;
    my @cigar;
    my $i = 0;
    while ($cigar_string =~ m/([0-9]+)([MIDNSHPX=])/g) {
	# skip hard-clipping
	die "ERROR: Unexpected CIGAR operation $2" if $2 eq "N" || $2 eq "P";
	if ($2 ne "H") {
	    push @cigar, [$1, $2];
	    $i += $1;
	}
    }
    return ($i, \@cigar);
}

sub min_max {
    my @array = sort {$a <=> $b} @_;
    return ($array[0], $array[$#array]);
}

sub reconstruct_alignment {
    my %alignment_object;
    my ($samfields_arrayref) = @_;
    ($alignment_object{cigar_total_length}, $alignment_object{cigar}) = &parse_cigar($samfields_arrayref->[5]);
    $alignment_object{readseq} = [split ("", uc($samfields_arrayref->[9]))];
    unshift @{$alignment_object{readseq}}, undef;
    $alignment_object{readqual} = [split ("", $samfields_arrayref->[10])];
    unshift @{$alignment_object{readqual}}, undef;
    $alignment_object{refseq} = "N" x $alignment_object{cigar_total_length};
    $alignment_object{refpos_start} = $samfields_arrayref->[3];
    $alignment_object{chrom} = $samfields_arrayref->[2];
    $alignment_object{direction} = $samfields_arrayref->[1] & 0x0010 ? "-" : "+";
    
    my @output;
    
    my @cigar_op; # 1-based, first element is not used
    my @cigar_refpos;
    my $cigar_pos_offset = 1;
    my @cigar_readseq;
    my @cigar_readqual;

    
    my $readbase_pos_offset = 0;
    my $refpos_pos_offset = $alignment_object{refpos_start};
    foreach my $cigarop (@{ $alignment_object{cigar} }) {
	my $len = $cigarop->[0];
	my $op = $cigarop->[1];
	map {$cigar_op[$_] = $op} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_readseq[$_] = $op eq "D" ? "*" : $alignment_object{readseq}->[$_ + $readbase_pos_offset]} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_readqual[$_] = $op eq "D" ? " " : $alignment_object{readqual}->[$_ + $readbase_pos_offset]} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_refpos[$_] = ($op eq "I" || $op eq "S") ? "*" : $refpos_pos_offset + $_ - $cigar_pos_offset} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	$readbase_pos_offset -= $len if $op eq "D";
	$refpos_pos_offset += $len if ($op ne "I" && $op ne "S");
	$cigar_pos_offset += $len;
    }

    ($alignment_object{left_alignpos}, $alignment_object{right_alignpos}) = &min_max(grep {defined $_ && $_ ne "*"} @cigar_refpos[1..$alignment_object{cigar_total_length}]);
    $alignment_object{cigar_op} = \@cigar_op;
    $alignment_object{cigar_refpos} = \@cigar_refpos;
    $alignment_object{cigar_readseq} = \@cigar_readseq;
    $alignment_object{cigar_readqual} = \@cigar_readqual;

    
    return \%alignment_object;
}

my %depthcache;
my $depthcache_count = 0;
sub depth {
    my ($region) = @_;
    return $depthcache{$region} if exists $depthcache{$region};
    my $output = `$samtools depth -d $max_samtools_depth -r $region $bam`;
    die "Error: unexpected output from samtools depth" unless $? == 0;
    my @outputlines = split("\n",$output);
    my @depth;
    foreach my $line (@outputlines) {
	my @fields = split (/\t/, $line);
	push @depth, $fields[2];
    }
    $depthcache{$region} = mean(@depth);
    $depthcache_count++;
    return $depthcache{$region};
}

sub mean {
    my $i = 0;
    my $sum = 0;
    map {$i++; $sum+=$_} @_;
    return $i == 0 ? -1 : $sum/$i
}