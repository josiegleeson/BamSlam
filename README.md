# BamSlam
This script was written for Oxford Nanopore Technologies' direct RNA sequencing data produced after mapping with minimap2 to the reference transcriptome. It will output important statistics and plots from your alignment data.

## Inputs:
To obtain a bam file align your FASTQ/FASTA files to the transcriptome with minimap2 (minimap2 -ax map-ont --sam-hit-only). You will also need a GTF annotation file. 

## Prerequisites:
<b>R packages:</b><br>
GenomicAlignments (from Bioconductor)<br>
GenomicFeatures (from Bioconductor)<br>
dplyr<br>
ggplot2<br>
data.table<br>
tidyr <br>
MASS <br>
viridis <br>

## Usage:
Download the Rscript from this repository and run it from the terminal as follows: <br>
Rscript BamSlam.R bam_file annotation_file out_prefix

## Outputs:
- A summary file of interesting metrics such as: percentage of full-length reads, median coverage fractions, percentage of reads that uniquely identify a transcript etc. <br>
- A full-length reads histogram. <br>
- 2D density plot of known transcript length vs coverage fractions. <br>
- A plot showing the number of reads with each number of secondary alignments. <br>
- A CSV file of all input alignments in case further analysis is required. <br>
- A CSV file summarised to median coverage of unique transcripts identified. <br>
