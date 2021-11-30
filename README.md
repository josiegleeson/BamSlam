# BamSlam
This script was written for Oxford Nanopore Technologies' direct RNA/cDNA sequencing data produced after mapping with minimap2 to the reference transcriptome. It will output a statistics file and plots from your alignment data. This script was used in: https://doi.org/10.1093/nar/gkab1129.

## Inputs:
To obtain a BAM file align your FASTQ/FASTA files to the transcriptome with minimap2 (minimap2 -ax map-ont --sam-hit-only). If minimap2 is not run with --sam-hit-only you should remove unmapped reads prior to running BamSlam to avoid slowing it down. You can also input a BAM file output from NanoCount: https://github.com/a-slide/NanoCount.

## Prerequisites:
<b>R packages:</b><br>
GenomicAlignments (from Bioconductor)<br>
dplyr<br>
tidyr<br>
tibble<br>
data.table<br>
ggplot2<br>
viridis <br>
hexbin<br>

## Usage:
Download/copy the Rscript from this repository and run it from the terminal as follows: <br>
Rscript BamSlam.R data_type bam_file out_prefix <br>
data_type: rna or cdna <br>

## Outputs:
- A summary CSV file of metrics such as: percentage of full-length reads, median coverage fractions, median read accuracy, median alignment length, number of transcripts identified, number of reads with no secondary alignments etc. <br>
- A full-length reads histogram (full-length cutoff/dashed line = 0.95). <br>
- A histogram density plot of known transcript length vs coverage fractions. <br>
- A bar chart showing number of secondary alignments. <br>
- A CSV file of input alignments (primary and secondary) in case further analysis is required. <br>
- A CSV file summarised to median coverage of unique transcripts identified. <br>
