ITDseek
=======
Detect _FLT3_ internal tandem duplication (_FLT3_ ITD) from amplicon sequencing reads (ITDseek) and simulate these reads (ITDsim)

[Download latest version in a ZIP package](https://github.com/tommyau/itdseek/zipball/master)

ITDseek for detection
---------------------
### Requirements, as tested on 64-bit CentOS 5.5
* SAMtools (version 0.1.19 tested)

### Usage
```bash
./itdseek.sh <sample.bam> <ref.fasta> <samtools> > sample.itdseek.vcf
```
- *sample.bam*: indexed BAM alignment file, generated from BWA-MEM with -M option (Mark shorter split hits as secondary). To index the BAM file: `samtools index sample.bam`
- *ref.fasta*: indexed hg19 reference genome in FASTA format. Only positions **chr13:28607161-28609590** will be considered, or `itdseek.sh` needs to be modified accordingly. To index the FASTA file: `samtools faidx ref.fasta`
- *samtools*: path to samtools executable. Use `samtools` if it is already included in a directory defined in `$PATH`
- *standard output (STDOUT)*: variant calls in [VCF version 4.1](http://samtools.github.io/hts-specs/VCFv4.1.pdf)

Example
```bash
# BWA-MEM with -M
bwa mem -R '@RG\tID:ITDsample\tSM:ITDsample' -M -t 12 ucsc.hg19.fasta ITDsample.R1.fastq ITDsample.R2.fastq | samtools view -bS - | samtools sort - ITDsample
# Index BAM
samtools index ITDsample.bam
# ITDseek
./itdseek.sh ITDsample.bam ucsc.hg19.fasta samtools > ITDsample.itdseek.vcf
```

### Interpretation of results
**Quality score** of called variants is defined as the **total number of sequencing reads with ITD detected** (forward and reverse combined).
The following is an example ITDseek output for a simulated 90bp ITD, with overall depth 2000X and VAF 50%. The quality score (sixth column `QUAL`) is 1000, corresponding to variant depth 1000X. Additional details are described in the last column (eighth column `INFO`), including forward and reverse ITD reads `DP2`, ITD length `LEN` and ITD sequence `SEQ`.
```bash
$ itdseek.sh ITD.100.90.bam ucsc.hg19.fasta samtools
##fileformat=VCFv4.1
##source=ITDseek
##reference=file://ucsc.hg19.fasta
##INFO=<ID=DP2,Number=2,Type=Integer,Description="# alt-foward and alt-reverse reads">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="length of ITD">
##INFO=<ID=SEQ,Number=1,Type=String,Description="sequence of ITD">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr13   28608300        .       C       CCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATC    1000    .       DP2=1000,0;LEN=90;SEQ=CATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATC
```
A simple filter based on quality score (i.e. number reads supporting ITD) is recommended for real experimental datasets, for example:
```bash
# Minimum quality score 50
./itdseek.sh ITDsample.bam ucsc.hg19.fasta samtools | awk '$1 ~ /^#/ || $6 >= 50' > ITDsample.itdseek.vcf
```


ITDsim for simulation
---------------------
### Requirements, as tested on 64-bit CentOS 5.5
* BioPerl (version 1.6.901 tested)

### Usage
```bash
./itdsim.pl | ./itdsim2fastq.pl
```
This generates the evaluation dataset described in the manuscript in the current directory, specifically 40401 pairs of FASTQ files (2x275bp at 2000X depth and ITD VAF 50%) representing 201 startings positions and 201 ITD lengths for the amplicon "FLT3.ITD.line.29.chr13.28607916.28608351_tile_2" in [Illumina TruSight Myeloid Panel](https://support.illumina.com/downloads/trusight-dna-amplicon-product-files.html).
Please refer to the comments in the source code of [itdsim.pl](https://github.com/tommyau/itdseek/blob/master/itdsim.pl) and [itdsim2fastq.pl](https://github.com/tommyau/itdseek/blob/master/itdsim2fastq.pl) for viewing/modifying simulation parameters.


Citation
--------
Au CH, Wa A, Ho DN, Chan TL and Ma ESK, 2016. [Clinical evaluation of panel testing by next-generation sequencing (NGS) for gene mutations in myeloid neoplasms.]((http://dx.doi.org/doi:10.1186/s13000-016-0456-8)) _Diagnostic Pathology_ 11:11 (doi:10.1186/s13000-016-0456-8)

License
-------
Source code released for non-commercial use only. For commercial use and other licensing enquiries, please contact Dr. Edmond S.K. Ma (<eskma@hksh.com>), Hong Kong Sanatorium and Hospital.