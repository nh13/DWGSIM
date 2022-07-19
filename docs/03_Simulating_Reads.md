# Simulating Reads

<!---toc start-->
  * [Overview](#overview)
  * [Error rates explained](#error-rates-explained)
  * [Read names explained](#read-names-explained)
  * [Mate pair or paired end modes](#mate-pair-or-paired-end-modes)
  * [Output mutations file](#output-mutations-file)
  * [Output FASTQ files](#output-fastq-files)

<!---toc end-->

## Overview

The `dwgsim` tool simulates reads from a reference genome FASTA.
It can be used to evaluate both mapping and variant calling.
See the `dwgsim_eval` tool for [Evaluating Mappings](04_Evaluating_Mappings.md).

Use `dwgsim -h` to see the full set of command line options (usage).

An example command is as follows:

```console
dwgsim -N 10000 -1 100 -2 100 -y 0 phix.fasta output
```

This will simulate 10,000 reads (`-N 10000`), which are paired end 2x100bp (`-1 100` for R1 and `-2 100` for R2), with no random reads (`-y 0`),
from the given genome FASTA (`phix.fasta`), producing output with prefix `output`.

The following output will be created:

| Name | Description |
| --- | --- |
| output.bfast.fastq.gz | Interleaved FASTQ containing both read one and read two |
| output.bwa.read1.fastq.gz | FASTQ containing only read one |
| output.bwa.read2.fastq.gz | FASTQ containing only read two |
| output.mutations.vcf | VCF containing simulated mutations |
| output.mutations.txt | TXT in a custom format containing simulated mutations (see [below](#output-mutations-file)) |


Notes:

- the longest supported insertion is 255bp
- The `-H` mode will simulate a haploid genome, whereas the default is to simulate a diploid genome. 

## Error rates explained 

The `-e` and `-E` options accept a uniform error rate (i.e. `-e 0.01` for 1%), or a uniformly increasing/decreasing error rate (i.e. `-e 0.01-0.1` for an error rate of 1% at the start of the read increasing to 10% at the end of the read).

## Read names explained

Read names are of the form:

```
 @<#1>_<#2>_<#3>_<#4>_<#5>_<#6>_<#7>_<#8>:<#9>:<#10>_<#11>:<#12>:<#13>_<#14>
```

| Field | Description |
| --- | --- |
| 1 | contig name (chromsome name) |
| 2 | start read 1 (one-based) |
| 3 | start read 2 (one-based) |
| 4 | strand read 1 (0 - forward, 1 - reverse) |
| 5 | strand read 2 (0 - forward, 1 - reverse) |
| 6 | random read 1 (0 - from the mutated reference, 1 - random) |
| 7 | random read 2 (0 - from the mutated reference, 1 - random) |
| 8 | number of sequencing errors read 1 (color errors for colorspace) |
| 9 | number of SNPs read 1 |
| 10 | number of indels read 1 |
| 11 | number of sequencing errors read 2 (color errors for colorspace) |
| 12 | number of SNPs read 2 |
| 13 | number of indels read 2 |
| 14 | read number (unique within a given contig/chromsome) |

Read 1 and read 2 correspond to the first and second reads from a paired-end/mate-pair read respectively.

## Mate pair or paired end modes

This utility can generate mate pair or paired end reads using the `-S` option.
By default, Illumina (nucleotide) data are paired end, and SOLiD (color space) data are mate pair.
For clarity, lets call the first end sequence E1 and the second end E2.

Paired end reads have the following orientation:

```
 5' E1 -----> ....             3'
 3'           .... <------- E2 5'
```

Above, the start co-ordinate of E1 is less than E2, with E1 and E2 reported on opposite strands.

Mate pair reads have following orientation

```
 5' E2 -----> .... E1 -------> 3'
 3'           ....             5'
```

Above, the start co-ordinate of E1 is greater than E2, with E1 and E2 reported on the same strand.

So for SOLiD mate pair reads, the R3 tag (E2) is listed before the F3 tag (E1).
For SOLiD paired end reads, the F3 tag (E1) is listed before the F5 tag (E2).

## Output mutations file

The locations of introduced mutations are given in a `<prefix>.mutations.txt` text file.
There are file columns:

1. the chromosome/contig name
2. the one-based position
3. the original reference base
4. the new reference base(s)
5. the variant strand(s)

SNPs are represented on one line, and in the case of heterozygous mutations, the new reference base is an [IUPAC](http://www.bioinformatics.org/sms/iupac.html) code.

```
 contig4   4   T   K   1
```

The above shows a heterozygous mutation at position 4 of contig4 on the first strand, mutating the T base to a heterozygous K (G or T) SNP.

Insertions are represented on one line, where the reference base is missing (indicated by a '-' in the third column).

```
 contig5   13   -   TAC   3
```

The above shows a homozygous insertion of TAC prior to position 13.

Each base of a deletion is represented on one line, where the new reference base is missing and represented by a '-'.

```
 contig6   22   A   -   2
```

The above shows a heterozygous deletion of T at position 22 on the second strand.
Multi-base deletions are show on consecutive lines.

```
 contig6   22   A   -   2
 contig6   23   C   -   2
```

The above shows a two base homozygous deletion of positions 22 and 23 on the second strand.

## Output FASTQ files

Three FASTQ files are produced, for use with BFAST (interleaved FASTQ) and BWA (one FASTQ per read end).

The FASTQ for BFAST is formatted so that the multi-end reads (paired end or mate pair) occur consecutively in the FASTQ (interleaved), with the read that is 5' of the other listed first.
For paired end reads, this means that E1 is always listed before E2.
For mate pair reads, this means that E2 is always listed before E1.

The FASTQs for BWA are split into two files, the first file for one end, the second file for the other, with the read that is 5' of the other in the first file.
For paired end reads, this means that E1 is in the first file and E2 is in the second file.
For mate pair reads, this means that E2 is in the first file and E1 is in the second file.

