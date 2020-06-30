# BamSlam
Filter your ONT BAM and get full-length transcripts.
This program has been designed for Oxford Nanopore Technologies' direct RNA sequencing data produced after mapping with minimap2 to the reference transcriptome. 

## Inputs:
You will need to obtain a bam file by aligning your FASTQ/FASTA files to the transcriptome with minimap2 (minimap2 -ax map-ont -N 100). You will also need a GTF annotation file. 

## Prerequisites:
<b>R packages:</b><br>
GenomicAlignments (from Bioconductor)<br>
GenomicFeatures (from Bioconductor)<br>
dplyr<br>
ggplot2<br>
data.table<br>
tidyr

## Usage:
Rscript BamSlam.R bam_file annotation_file out_prefix

## Outputs:
- A summary file of the filtered alignment metrics such as: number and percentage of full-length reads, median coverage fractions etc. <br>
- A full-length reads histogram. <br>
- 2D density plot of known transcript length vs coverage fractions. <br>
- A CSV file of all input alignments in case further analysis is required. <br>
