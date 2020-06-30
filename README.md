# BamSlam
This script was written for Oxford Nanopore Technologies' direct RNA sequencing data produced after mapping with minimap2 to the reference transcriptome. It will output important statistics ad output plots from your alignment data.

## Inputs:
To obtain a bam file align your FASTQ/FASTA files to the transcriptome with minimap2 (minimap2 -ax map-ont -N 100). Then discard unampped and supplementary alignments with samtools (samtools view -hb -F 2052 in.bam > out.bam). You will also need a GTF annotation file. 

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
- A CSV file of all input alignments in case further analysis is required. <br>
