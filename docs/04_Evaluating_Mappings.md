# Evaluating Mappings

<!---toc start-->
  * [Overview](#overview)
  * [Key Options](#key-options)
  * [Output](#output)

<!---toc end-->

## Overview

The `dwgsim_eval` tool evaluates the mappings from reads produced by `dwgsim`.
See the `dwgsim` tool for [Simulating Reads](03_Simulating_Reads.md).

Use `dwgsim_eval -h` to see the full set of command line options (usage).

An example command is as follows:

```console
dwgsim_eval mapped.bam
```

A read is correctly mapped if the mapped start position is within `INT` bases (`-g`) of the simulated read's start position.
Otherwise, it is counted as incorrectly mapped.
Furthermore, random reads created by `dwgsim` are expected to be unmapped (counted as false mappings otherwise).

## Key Options

- Use the `-c` option if the reads are from SOLID
- Use the `-b` option if the alignments are from BWA when mapping SOLiD reads
- Use the `-z` option if the input reads are expected to be single end (not paired end)
- Use the `-s` option if the input is a SAM file (default: BAM)
- Use the `-n` option to specify the number of reads if only a subset of the input reads were used (affects the denominator)


## Output

The `dwgsim_eval` tool outputs a table with the following columns:

| Column | Description |
| --- | --- |
| `thr` | The alignment threshold (see `-a`) |
| `mc` | The number of reads mapped correctly that should be mapped at the threshold |
| `mi` | The number of reads mapped incorrectly that should be mapped be mapped at the threshold |
| `mu` | The number of reads unmapped that should be mapped be mapped at the threshold |
| `um` | The number of reads mapped that should be unmapped be mapped at the threshold |
| `uu` | The number of reads unmapped that should be unmapped be mapped at the threshold |
| `mc + mi + mu + um + uu` | The total number of reads that should be unmapped be mapped at the threshold |
| `mc'` | The number of reads mapped correctly that should be mapped at or greater than the threshold |
| `mi'` | The number of reads mapped incorrectly that should be mapped be mapped at or greater than the threshold |
| `mu'` | The number of reads unmapped that should be mapped be mapped at or greater than the threshold |
| `um'` | The number of reads mapped that should be unmapped be mapped at or greater than the threshold |
| `uu'` | The number of reads unmapped that should be unmapped be mapped at or greater than the threshold |
| `mc' + mi' + mu' + um' + uu'` | The total number of reads that should be unmapped be mapped at or greater than the threshold |
| `mc / (mc + mi + mu)` | Sensitivity at the threshold. I.e. the fraction of reads that should be mapped that are mapped correctly. |
| `mc / (mc + mi)` | Positive predictive value at the threshold. I.e. The fraction of mapped reads that are mapped correctly. |
| `um / (um + uu)` | False discovery rate at the threshold.  I.e. The fraction of random reads that are mapped. |
| `mc' / (mc' + mi' + mu')` | Sensitivity at or greater than the threshold. I.e. the fraction of reads that should be mapped that are mapped correctly. |
| `mc' / (mc' + mi')` | Positive predictive value at or greater than the threshold. I.e. The fraction of mapped reads that are mapped correctly. |
| `um' / (um' + uu')` | False discovery rate at or greater than the threshold.  I.e. The fraction of random reads that are mapped. |

"At or greater than the threshold" tells us what our sensitivity, PPV, and FDR would be if we filtered based on that threshold.
