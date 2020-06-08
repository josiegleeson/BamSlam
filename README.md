# BamSlam
Filter and extract data from ONT minimap2 RNASeq alignments.

## Inputs:
<br>
You will need to align your FASTQ/FASTA files to the transcriptome with minimap2 (minimap2 -ax map-ont -N 100). Then create a BAM file of all primary and secondary alignments. You will also need a GTF/GFF annotation file. 
<br>
<br>
## Prerequisites:
<br>

<br>
<br>
## Installation:
<br>

<br>
<br>
## Usage:
<br>
bash BamSlam bam_file annotation_file out_prefix
<br>
<br>
## Outputs:
<br>
- A BAM file of filtered primary and secondary alignments. You can use this as input into Salmon for quantification.
- A CSV file of all input alignments in case further analysis is required.
- A summary file filtered alignment metrics such as: number and percentage of full-length reads, median coverage fractions etc. 
