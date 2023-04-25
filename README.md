# BamSlam
This script was written for Oxford Nanopore Technologies long-read direct RNA/cDNA sequencing data produced after mapping with minimap2 to the reference transcriptome. It will output a summary file and plots from your alignment data. This script was used in: https://doi.org/10.1093/nar/gkab1129.

## Inputs:
To obtain a BAM file align your FASTQ/FASTA files to the transcriptome with minimap2.
```
minimap2 -ax map-ont --sam-hit-only transcriptome.fasta input.fastq > alignments.sam
```
If minimap2 is not run with --sam-hit-only you should remove unmapped reads prior to running BamSlam to avoid slowing it down. You can also input a BAM file output from NanoCount: https://github.com/a-slide/NanoCount.

## Prerequisites:
<b>R packages:</b>
- GenomicAlignments (Bioconductor)
- dplyr
- tidyr
- tibble
- data.table
- ggplot2
- viridis
- hexbin

## Usage:
Download/copy the Rscript from this repository and run it from the terminal as follows:

```
Rscript BamSlam.R [DATA_TYPE] [BAM_FILE] [OUT_PREFIX]
Rscript BamSlam.R rna undiff1_5Y.bam undiff1

DATA_TYPE, enter either: cdna, rna
BAM_FILE, a BAM file of alignments to the transcriptome
OUT_PREF, output file prefix
```

## Outputs:
- A summary CSV file of metrics such as: percentage of full-length reads, median coverage fractions, median read accuracy, median alignment length, number of transcripts identified, number of reads with no secondary alignments etc.
- A full-length reads histogram (full-length cutoff/dashed line = 0.95). 
- A histogram density plot of known transcript length vs coverage fractions. 
- A bar chart showing number of secondary alignments.
- A CSV file of input alignments (primary and secondary) in case further analysis is required.
- A CSV file summarised to median coverage of unique transcripts identified.
