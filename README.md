ITDseek
=======
Detect _FLT3_ internal tandem duplication (_FLT3_ ITD) from amplicon sequencing reads (ITDseek) and simulate these reads (ITDsim)

[Download latest version in a ZIP package](https://github.com/tommyau/itdseek/zipball/master)

_FLT3_ ITD detection performance of **ITDseek, Pindel, GATK HaplotypeCaller and Samtools** based on 40,401 simulated samples of ITD length up to 201bp at different positions of an amplicon (Au _et al. Diagnostic Pathology_ 2016 11:11 [Figure 2](http://link.springer.com/article/10.1186/s13000-016-0456-8/fulltext.html#Fig2))

![FLT3 ITD detection performance of ITDseek, Pindel, GATK HaplotypeCaller and Samtools](http://static-content.springer.com/image/art%3A10.1186%2Fs13000-016-0456-8/MediaObjects/13000_2016_456_Fig2_HTML.gif)

ITDseek for detection
---------------------
### Requirements, as tested on 64-bit CentOS 5.5
* [SAMtools](http://www.htslib.org/download/) version at least 1.3 ([to support depth >8000x](https://github.com/samtools/samtools/pull/322))

### Usage
```bash
./itdseek.sh <sample.bam> <ref.fasta> <samtools> > sample.itdseek.vcf
```
- *sample.bam*: indexed BAM alignment file, generated from BWA-MEM with -M option (Mark shorter split hits as secondary). To index the BAM file: `samtools index sample.bam`
- *ref.fasta*: indexed hg19 reference genome in FASTA format. Only positions **chr13:28607161-28609590** will be considered, or `itdseek.sh` needs to be modified accordingly. To index the FASTA file: `samtools faidx ref.fasta`
- *samtools*: path to samtools executable. Use `samtools` if it is already included in a directory defined in `$PATH`
- *standard output (STDOUT)*: variant calls in [VCF version 4.1](http://samtools.github.io/hts-specs/VCFv4.1.pdf)

Examples
```bash
# Provided example data
# 100bp ITD leading to soft-clipping in BAM file
./itdseek.sh examples/ITD.100.100.bam ucsc.hg19.fasta samtools > ITD.100.100.itdseek.vcf
# 20bp ITD leading to insertion in BAM file
./itdseek.sh examples/ITD.100.20.bam ucsc.hg19.fasta samtools > ITD.100.20.itdseek.vcf

# Your own data
# BWA-MEM with -M
bwa mem -R '@RG\tID:ITDsample\tSM:ITDsample' -M ucsc.hg19.fasta ITDsample.R1.fastq ITDsample.R2.fastq | samtools sort > ITDsample.bam
# Index BAM
samtools index ITDsample.bam
# ITDseek
./itdseek.sh ITDsample.bam ucsc.hg19.fasta samtools > ITDsample.itdseek.vcf
```

### Interpretation of results
**Quality score** of called variants is defined as the **total number of sequencing reads with ITD detected** (forward and reverse combined).
The following is ITDseek output for a simulated 100bp ITD dataset included in the folder examples/, with overall depth 2000X and VAF 50%. The quality score (sixth column `QUAL`) is 1000, corresponding to variant depth 1000X. Additional details are described in the last column (eighth column `INFO`), including forward and reverse ITD reads `DP2`, ITD length `LEN` and ITD sequence `SEQ`.
```bash
$ ./itdseek.sh examples/ITD.100.100.bam ucsc.hg19.fasta samtools
##fileformat=VCFv4.1
##source=ITDseekV1.2
##reference=file://ucsc.hg19.fasta
##INFO=<ID=DP2,Number=2,Type=Integer,Description="# alt-foward and alt-reverse reads">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="length of ITD">
##INFO=<ID=SEQ,Number=1,Type=String,Description="sequence of ITD">
##INFO=<ID=CLIPPING,Number=0,Type=Flag,Description="ITD is detected as soft-clipping">
##INFO=<ID=INSERTION,Number=0,Type=Flag,Description="ITD is detected as insertion">
##INFO=<ID=VAF,Number=1,Type=Float,Description="ITD allele fraction">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr13	28608310	.	G	GATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCG	2000	.	DP2=1000,1000;LEN=100;SEQ=ATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCG;CLIPPING;VAF=0.50
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
Au CH, Wa A, Ho DN, Chan TL and Ma ESK, 2016. [Clinical evaluation of panel testing by next-generation sequencing (NGS) for gene mutations in myeloid neoplasms.](http://link.springer.com/article/10.1186/s13000-016-0456-8/fulltext.html) _Diagnostic Pathology_ 11:11 (doi:10.1186/s13000-016-0456-8)

License
-------
Source code released for non-commercial use only. For commercial use and other licensing enquiries, please contact Dr. Edmond S.K. Ma (<eskma@hksh.com>), Hong Kong Sanatorium and Hospital.