![](http://nugendata.com/images/nugen_logo_noedge.png)  
[www.nugen.com](http://www.nugen.com)

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)

System Requirements
-----------------------------
- samtools required [samtools page](https://github.com/samtools/samtools/releases/) Tested on 0.1.19, 1.2.0, 1.3.0
- python2.7 [anaconda page](http://continuum.io/downloads) Tested on 2.7.7 
- GNU-coreutils gzip,sed,cut,grep etc. Standard on Ubuntu 12.04+.
- Aligner: STAR or BWA. (Tophat and Tophat 2 do not set SAM flags used by NuDUP and is not supported.) If a spliced aligner is not required, bowtie/bowtie2 are compatibile with `nudup`.

Usage
-----------------------------
Run `python nudup.py -h` for latest usage.

```
usage: nudup.py [-2] [-f INDEX.fq|READ.fq] [-o OUT_PREFIX] [-s START]
                [-l LENGTH] [-v] [-h]
                IN.sam|IN.bam

Marks/removes PCR introduced duplicate molecules based on the molecular tagging
technology used in NuGEN products.

For SINGLE END reads, duplicates are marked if they fulfill the following
criteria: a) start at the same genomic coordinate b) have the same strand
orientation c) have the same molecular tag sequence. The read with the
highest mapping quality is kept as the non-duplicate read.

For PAIRED END reads, duplicates are marked if they fulfill the following
criteria: a) start at the same genomic coordinate b) have the same template
length c) have the same molecular tag sequence. The read pair with the highest
mapping quality is kept as the non-duplicate read.

Here are the two cases for running this tool:

- CASE 1 (Standard): User supplies two input files,
 1) SAM/BAM file that a) ends with .sam or .bam extension b) contains unique
    alignments only
 2) FASTQ file (ie: the Index FASTQ) that contains the molecular tag sequence
    for each read name in the corresponding SAM/BAM file as either a) the read
    sequence or b) in the FASTQ entry header name. If the index FASTQ read
    length is 6, 8, 12, 14, or 16nt long as expected for NuGEN products, the
    molecular tag sequence to be extracted from the read according to -s and -l
    parameters, otherwise the molecular tag will be extracted from the header
    of the FASTQ entry.

- CASE 2 (Runtime Optimized): User supplies only one input file,
 1) SAM/BAM file that a) ends with .sam or .bam extension b) contains unique
    alignments only c) is sorted d) has a fixed length sequence containing the
    molecular tag appended to each read name.

Author: Anand Patel
Contact: NuGEN Technologies Inc., techserv@nugen.com

Input:
  IN.sam|IN.bam         input sorted/unsorted SAM/BAM containing only unique
                        alignments (sorted required for case 2 detailed above)

Options:
  -2, --paired-end      use paired end deduping with template. SAM/BAM
                        alignment must contain paired end reads. Degenerate
                        read pairs (alignments for one read of pair) will be
                        discarded.
  -f INDEX.fq|READ.fq   FASTQ file containing the molecular tag sequence for
                        each read name in the corresponding SAM/BAM file
                        (required only for CASE 1 detailed above)
  -o OUT_PREFIX, --out OUT_PREFIX
                        prefix of output file paths for sorted BAMs (default
                        will create prefix.sorted.markdup.bam,
                        prefix.sorted.dedup.bam, prefix_dup_log.txt)
  -s START, --start START
                        position in index read where molecular tag sequence
                        begins. This should be a 1-based value that counts in
                        from the 3' END of the read. (default = 6)
  -l LENGTH, --length LENGTH
                        length of molecular tag sequence (default = 6)
  -T TEMP_DIR           directory for reading and writing to temporary files
                        and named pipes (default: /tmp)
  --old-samtools        required for compatibility with samtools sort style in
                        samtools versions <=0.1.19
  --rmdup-only          required for only outputting duplicates removed file
  -v, --version         show program's version number and exit
  -h, --help            show this help message and exit
```

Support
-----------------------------  
For questions, contact NuGEN Technologies Technical Support techserv@nugen.com
