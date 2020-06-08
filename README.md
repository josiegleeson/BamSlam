# BamSlam
Filter your ONT BAM and get full-length transcripts.

## What is BamSlam?
This program has been designed for Oxford Nanopore Technologies' direct RNA sequencing data produced after mapping with minimap2 to the reference transcriptome. In many cases, the primary alignment of a read may be incorrect so simply using primary alignments for downstream analysis may give misleading results (ref: Soneson et al. 2019). When a read has a primary and secondary alignment with the same chaining score, minima2p2 assings the primary alignment randomly. In the resulting alignment file, there are many reads with multiple secondary alignments with identical alignment scores to the primary. BamSlam was created to filter secondary alignments to retain only those which are highly similar to the primary alignment.

## How does it work?
Due to the mechanism of dRNA sequeuncing, reads always enter the nanopore at the 3 prime end first. So we know the 3 prime ends are the true beginning of the reads. The first step of BamSlam filters the alignments to keep only those that have a 3' site within 50nt of the transcript 3' end. Next, different parameters for every alignment are compared. Each alignment within a read is compared, keeping only alignments with an alignment score that is at least 98% as good as the best alignment score for that read. This is then followed for alignments with an aligned fraction and accuracy of at least 90% of the best alignment.

## Inputs:
You will need to align your FASTQ/FASTA files to the transcriptome with minimap2 (minimap2 -ax map-ont -N 100). Then create a BAM file of all primary and secondary alignments. You will also need a GTF/GFF annotation file. 

## Installation:
Download the scripts by cloning this repository. The scripts can then be run directly given dependencies are already installed. 

## Prerequisites:
<b>Samtools</b> (tested with  v1.10)<br>
<b>R</b> (tested with v3.6.1)<br>
<b>R packages:</b><br>
GenomicAlignments (from Bioconductor)<br>
GenomicFeatures (from Bioconductor)<br>
dplyr<br>
ggplot2<br>
data.table<br>
tidyr

## Usage:
BamSlam bam_file annotation_file out_prefix

## Outputs:
- A BAM file of filtered primary and secondary alignments. You can use this as input into Salmon for quantification. <br>
- A CSV file of all input alignments in case further analysis is required. <br>
- A summary file of the filtered alignment metrics such as: number and percentage of full-length reads, median coverage fractions etc. <br>
