# BamSlam
Filter and extract data from ONT minimap2 RNASeq alignments.

## Inputs:
You will need to align your FASTQ/FASTA files to the transcriptome with minimap2 (minimap2 -ax map-ont -N 100). Then create a BAM file of all primary and secondary alignments. You will also need a GTF/GFF annotation file. 
<br>
## Prerequisites:
<b>Samtools<br>
R and the following libraries:</b><br>
GenomicAlignments (from Bioconductor)<br>
GenomicFeatures (from Bioconductor)<br>
dplyr<br>
ggplot2<br>
data.table<br>
tidyr

## Installation:


## Usage:
BamSlam bam_file annotation_file out_prefix

## Outputs:
- A BAM file of filtered primary and secondary alignments. You can use this as input into Salmon for quantification. <br>
- A CSV file of all input alignments in case further analysis is required. <br>
- A summary file of the filtered alignment metrics such as: number and percentage of full-length reads, median coverage fractions etc. <br>
