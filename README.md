nudup -- NuGEN's tool for deduplicating reads using N6 sequence information
=============================

Software used
-----------------------------
- samtools-0.1.18 or higher [samtools page](http://samtools.sourceforge.net/) Tested on 0.1.18
- python2.7 [anaconda page](http://continuum.io/downloads) Tested on 2.7.7 
- GNU-coreutils gzip,sed,cut,grep etc. Standard on Ubuntu 12.04+.

Usage
-----------------------------
Run `python nudup.py -h` for latest usage.

usage: nudup.py &#91;-2&#93; &#91;-f INDEX.fq|READ.fq&#93; &#91;-o OUT&#95;PREFIX&#93; &#91;-s START&#93;
                &#91;-l LENGTH&#93; &#91;-v&#93; &#91;-h&#93;
                IN.sam|IN.bam

Marks/removes PCR duplicates using the N6 sequence technology used in
NuGEN's Ovation Target Enrichment libraries.

For single end reads, duplicates are marked if they fulfill the following
criteria: a) start in the same location b) have the same orientation/strand
c) have the same N6 sequence. The read with the highest mapping quality is kept
as the non-duplicate read. For paired end reads, duplicates are marked if they
fulfill the following criteria: a) start in the same location b) have the same
template length c) have the same N6 sequence. The read pair with the highest
mapping quality is kept as the non-duplicate read.

The runtime for marking and removal of duplicates depends on the format of the
inputs. Here are the two cases for running this tool:

- Case 1 (faster runtime): User supplies pre-built SAM/BAM alignment file that
  is sorted and has fixed length N6 sequence appended to the read name
- Case 2 (slower runtime): User supplies SAM/BAM alignment file and FASTQ file
  containing N6 sequence. SAM/BAM records do not have a N6 sequence, but the N6
  sequence is a substring in the read or title of the FASTQ file. FASTQ read
  names must be unique, and each SAM/BAM record must have corresponding read
  name in the FASTQ file.

Author: Anand Patel
Contact: NuGEN Technologies Inc., techserv@nugen.com

Input:
  IN.sam|IN.bam         input sorted/unsorted SAM/BAM

Options:
  -2, --paired-end      use paired end deduping with template. SAM/BAM
                        alignment must contain paired end reads.
  -f INDEX.fq|READ.fq   FASTQ file paired with input SAM/BAM (REQUIRED if
                        SAM/BAM does not have N6 sequence appended to
                        qname/read title)
  -o OUT&#95;PREFIX, --out OUT&#95;PREFIX
                        prefix of output file paths for sorted BAMs (default
                        will create prefix.sorted.markdup.bam,
                        prefix.sorted.dedup.bam, prefix_dup_log.txt)
  -s START, --start START
                        position in index read where N6 sequence begins. This
                        should be a 1-based value that counts in from the END
                        of the read. (default = 6)
  -l LENGTH, --length LENGTH
                        length of N6 sequence (default = 6)
  -v, --version         show program's version number and exit
  -h, --help            show this help message and exit

