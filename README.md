# BamSlam
This script was written for long-read Oxford Nanopore Technologies direct RNA/cDNA sequencing data. It uses BAM files produced after mapping with minimap2 to the reference transcriptome. It will output a summary file and plots from your aligned reads. This script was used in: https://doi.org/10.1093/nar/gkab1129.

## Inputs:
To obtain a BAM file align your FASTQ/FASTA files to the transcriptome with minimap2.
```
minimap2 -ax map-ont --sam-hit-only transcriptome.fasta reads.fastq | samtools view -bh > aligned_reads.bam
```
If minimap2 is not run with --sam-hit-only you should remove unmapped reads prior to running BamSlam to avoid slowing it down. You can also input a BAM file output from NanoCount: https://github.com/a-slide/NanoCount.

## Requirements:
R (tested with v4.2.2)
R packages:
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
The script takes approximately 5 minutes per million reads.

## Outputs:
- A summary CSV file of your alignments.
- A CSV file of all input alignments (primary and secondary) for downstream analysis.
- A CSV file of each unique transcript identified in the data with its corresponding median read coverage.
- A histogram of read coverages (full-length reads cutoff/dashed line = 0.95).
- A histogram of known transcript lengths and the number of reads assigned. 
- A histogram density plot of known transcript length vs coverage fractions. 
- A bar chart of the secondary alignments.
- Density plot of the read accuracies.

# Plots
<img src="https://github.com/josiegleeson/BamSlam/assets/30969357/c1d2d5f6-066f-407a-b686-8028e3bd06f3" width="100" height="100">
![new_coverage_fraction](https://github.com/josiegleeson/BamSlam/assets/30969357/c1d2d5f6-066f-407a-b686-8028e3bd06f3)
![new_transcript_length_distribution](https://github.com/josiegleeson/BamSlam/assets/30969357/a785bcb4-4472-41e8-ad1c-04481ba085f0)
![new_density](https://github.com/josiegleeson/BamSlam/assets/30969357/9274dc53-dc61-4cef-8158-cbe5ac89984b)
![new_sec_alns](https://github.com/josiegleeson/BamSlam/assets/30969357/2634921a-d755-48da-91ec-39966e0de4a6)
![new_accuracy](https://github.com/josiegleeson/BamSlam/assets/30969357/5047fa9f-ac44-45a5-9c29-8fccfffad3dd)


## Detailed description of summary file:
- Total number of reads
- Number of reads representing full-length transcripts (reads with coverage fractions > 0.95)
- Percentage of reads representing full-length transcripts
- Median read coverage fraction (primary alignments)
- Median alignment length (primary alignments)
- Median accuracy (primary alignments)
- Number of reads with no secondary alignments
- Percentage of reads with no secondary alignments
- Total number of distinct transcripts identified in the data
- Median coverage fraction of all unique transcripts (median value of the median read coverages per transcript)
- Median length of all unique transcripts identified (median length in nt of the annotated length of transcripts identified)


