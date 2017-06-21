#!/bin/bash
# Chun Hang AU (chau@hksh.com) (https://github.com/tommyau/itdseek)
# ITDseek - FLT3 ITD detection in amplicon sequencing reads
# Example usage:
# itdseek.sh ITD.bam ucsc.hg19.fasta samtools

SCRIPT_PATH=$(readlink -f $0)
SCRIPT_DIR=$(dirname $SCRIPT_PATH)
SAMTOOLS_VERSION_REQUIRED=1.3
# argument check
if [[ "$#" -eq 3 && "$1" != "" && "$2" != "" && "$3" != "" && -f "$1" && -f "$2" ]]; then
    BAM="$1"
    REFSEQ="$2"
    SAMTOOLS="$3"
    BAMbn=$(basename "$BAM")
else
    echo $(date) "[ERROR] Please define arguments as: itdseek.sh [bam] [refseq] [samtools]. Example: ./itdseek.sh ITD.bam ucsc.hg19.fasta samtools"
    exit 1
fi

# check samtools version
function version { echo "$@" | cut -f1 -d"+" | awk -F. '{ printf("%03d%03d%03d\n", $1,$2,$3); }'; }
"$SAMTOOLS" --version-only >/dev/null 2>&1 || { echo >&2 "[ERROR] SAMtools (provided path: $SAMTOOLS) is not running properly. Please check the path and/or version (at least $SAMTOOLS_VERSION_REQUIRED)"; exit 1; }
SAMTOOLS_VERSION=`"$SAMTOOLS" --version-only`
if [ "$(version "$SAMTOOLS_VERSION")" -lt "$(version "$SAMTOOLS_VERSION_REQUIRED")" ]; then
     echo >&2 "[ERROR] SAMtools version ($SAMTOOLS_VERSION) is not supported (supported version: at least $SAMTOOLS_VERSION_REQUIRED)."
     exit 1
fi

# pipeline SAM alignments to itdseek.pl
"$SAMTOOLS" view "$BAM" chr13:28607161-28609590 | "${SCRIPT_DIR}/itdseek.pl" --refseq "$REFSEQ" --samtools "$SAMTOOLS" --bam "$BAM"
