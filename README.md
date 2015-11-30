ITDseek
=======
Detect FLT3 internal tandem duplication (FLT3 ITD) from amplicon sequencing reads (ITDseek) and simulate these reads (ITDsim)

[Download latest version in a ZIP package](https://github.com/tommyau/itdseek/zipball/master)

ITDseek for detection
---------------------
### Requirements, as tested on 64-bit CentOS 5.5
* SAMtools (version 0.1.19 tested)

### Usage
```bash
./itdseek.sh <sample.bam> <ref.fasta> <samtools> > sample.itdseek.vcf
```
- *<sample.bam>*: indexed BAM alignment file, generated from BWA-MEM with -M option (Mark shorter split hits as secondary). To index the BAM file: `samtools index sample.bam`
- *<ref.fasta>*: indexed hg19 reference genome in FASTA format. Only positions **chr13:28607161-28609590** will be considered, or `itdseek.sh` needs to be modified accordingly. To index the FASTA file: `samtools faidx ref.fasta`
- *<samtools>*: path to samtools executable. Use `samtools` if it is already included in a directory defined in `$PATH`
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


ITDsim for simulation
---------------------
### Requirements, as tested on 64-bit CentOS 5.5
* BioPerl (version 1.6.901 tested)

### Usage
```bash
./itdsim.pl | ./itdsim2fastq.pl
```
This generates the evaluation dataset described in the manuscript in the current directory, specifically 40401 pairs of FASTQ files (2x275bp at 2000X depth and ITD VAF 50%) representing 201 startings positions and 201 ITD lengths for the amplicon "FLT3.ITD.line.29.chr13.28607916.28608351_tile_2" in [Illumina TruSight Myeloid Panel](https://support.illumina.com/downloads/trusight-dna-amplicon-product-files.html)


Citation
--------
Au CH, Wa A, Ho DN, Chan TL and Ma ESK, 2015. Clinical evaluation of panel testing by next-generation sequencing for gene mutations in myeloid neoplasms *(submitted)*

License
-------
Source code released for non-commercial use only. For commercial use and other licensing enquiries, please contact Dr. Edmond S.K. Ma (<eskma@hksh.com>), Hong Kong Sanatorium and Hospital.