ITDseek
=======
Detect FLT3 internal tandem duplication (FLT3 ITD) from amplicon sequencing reads (ITDseek) and simulate these reads (ITDsim)

ITDseek for detection
---------------------
### Requirements, as tested on 64-bit CentOS 5.5
* SAMtools (version 0.1.19 tested)

### Procedure
```
./itdseek.sh [bam] [refseq] [samtools]
```

Example
```
./itdseek.sh ITD.bam ucsc.hg19.fasta samtools
```


ITDsim for simulation
---------------------
### Requirements, as tested on 64-bit CentOS 5.5
* BioPerl (version 1.6.901 tested)

### Procedure
```
./itdsim.pl | ./itdsim2fastq.pl
```
