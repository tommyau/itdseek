#!/bin/bash
# Chun Hang AU (chau@hksh.com) (https://github.com/tommyau/itdseek)
# ITDseek - FLT3 ITD detection in amplicon sequencing reads
# Example usage:
# itdseek.sh ITD.bam ucsc.hg19.fasta samtools

SCRIPT_PATH=$(readlink -f $0)
SCRIPT_DIR=$(dirname $SCRIPT_PATH)
# argument check
if [[ "$#" -eq 3 && "$1" != "" && "$2" != "" && "$3" != "" ]]; then
    BAM="$1"
    REFSEQ="$2"
    SAMTOOLS="$3"
else
    echo $(date) "[ERROR] Please define arguments as: itdseek.sh [bam] [refseq] [samtools]. Example: itdseek.sh ITD.bam ucsc.hg19.fasta samtools"
    exit 1
fi

# pipeline SAM alignments to itdseek.pl
$SAMTOOLS view $BAM chr13:28607161-28609590 | ${SCRIPT_DIR}/itdseek.pl --refseq $REFSEQ --samtools $SAMTOOLS