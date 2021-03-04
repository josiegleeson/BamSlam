# Replace sample1 with your file name or sample name
# Usage: Rscript BamSlam.R rna/cdna yourfile.bam ouputprefix

# GenomicAlignments/Features package is from bioconductor, need to install bioconductor then run:
# BiocManager::install("GenomicFeatures")

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  type <- args[1]
  bamfile <- args[2]
  output <- args[3]
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(dplyr, warn.conflicts = FALSE)
    library(tibble)
    library(tidyr)
    library(ggplot2)
    library(viridis)
  })
    
  options(dplyr.summarise.inform = FALSE)
  
  # Import bam
  bam <- readGAlignments(bamfile, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                              what = c("qname","flag","mapq")))
  message("Imported bam file")
  
  # Expand CIGAR strings
  cigaropts = cigarOpTable(bam@cigar)
  
  for (col in colnames(cigaropts)) {
    mcols(bam)[[paste0("nbr", col)]] = cigaropts[, col]
  }
  
  message("Expanded CIGAR strings")
  
  bam_1 <- as.data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::select(-cigar) %>% 
    dplyr::select(-njunc)
  
  if (type == "cdna") {
    bam_data <- subset(bam_1, flag == 0 | flag == 16 | flag == 256 | flag == 272)
  } else if (type == "rna") {
    bam_data <- subset(bam_1, flag == 0 | flag == 256)
  } else {
    print("Sequencing type missing or invalied. Please enter either: cdna rna")
  }
  
  bam_data <- bam_data %>% 
    dplyr::mutate(alignedLength = nbrM + nbrI) %>% 
    dplyr::mutate(readLength = nbrS + nbrH + nbrM + nbrI) %>% 
    dplyr::mutate(alignedFraction = alignedLength/readLength) %>% 
    dplyr::mutate(accuracy=(nbrM+nbrI+nbrD-NM)/(nbrM+nbrI+nbrD))
  
  lengths <- as.data.frame(bam@seqinfo) %>% 
    dplyr::select(-isCircular) %>% 
    dplyr::select(-genome) %>% 
    tibble::rownames_to_column("reference")
  
  # Merge known lengths and calculate transcript coverage
  bam_data <- merge(bam_data, lengths, by.x="seqnames", by.y="reference", all.x=TRUE)
  
  bam_data <- bam_data %>% 
    dplyr::mutate(coverage=width/seqlengths)
  
  bam_data$coverage <- as.numeric(bam_data$coverage)
  
  bam_data <- bam_data %>% 
    tidyr::drop_na(coverage)      
  
  message("Calculated transcript coverages") 

  alignments <- bam_data %>% 
    group_by(qname) %>% 
    summarise(nbrSecondary = n()-1) %>% 
    rename(read = qname)
  
  bam_data <- merge(bam_data, alignments, by.x="qname", by.y="read", all.x=TRUE)
  
  alignments <- alignments %>% 
    group_by(nbrSecondary) %>% 
    summarise(total = n()) %>%
    mutate(prop = total / sum(total))
    
  bam_data$nbrSecondary <- as.factor(bam_data$nbrSecondary)
  alignments$nbrSecondary <- as.factor(alignments$nbrSecondary)
  
  # Export files
  bam_export <- subset(bam_data, select=c("qname", "seqnames", "flag", "mapq", "AS", "alignedLength", "readLength", "alignedFraction", "accuracy", "seqlengths", "coverage", "nbrSecondary"))
  write.csv(bam_export, file = paste0(output, "_data.csv"), sep=",", quote=F, col.names = T, row.names=F) 
  
  bam_primary <- bam_data %>% 
    subset(flag == 0 | flag == 16)
    
  bam_per_unique_transcript <- bam_primary %>% 
    dplyr::group_by(seqnames) %>% 
    summarise(coverage = median(coverage, na.rm = TRUE))
  
  length_per_unique_transcript <- bam_primary %>% 
    dplyr::group_by(seqnames) %>% 
    summarise(seqlengths = max(seqlengths))
  
  write.csv(bam_per_unique_transcript, file = paste0(output, "_transcript_level_data.csv"), sep=",", quote=F, col.names = T, row.names=F)  
  
  a <- sum(bam_primary$coverage > 0.95)
  b <- nrow(bam_primary)
  c <- a/nrow(bam_primary)*100
  d <- median(bam_primary$coverage)
  e <- median(bam_primary$accuracy)*100
  f <- sum(bam_data$nbrSecondary == 0)
  g <- f/nrow(bam_primary)*100
  h <- nrow(bam_per_unique_transcript)
  i <- median(bam_per_unique_transcript$coverage)
  j <- median(length_per_unique_transcript$seqlengths)
  
  # Make stats
  metric <- c(
    "Number of reads representing full-length transcripts:",
    "Out of total primary alignments:",
    "Percentage of reads representing full-length transcripts:",
    "Median coverage fraction of transcripts (primary alignments):",
    "Median accuracy of primary alignments:",
    "Number of reads with no secondary alignments:",
    "Percentage of reads with no secondary alignments:",
    "Number of unique transcripts identified:",
    "Median coverage fraction of all unique transcripts:",
    "Median length of all unique transcripts identified:") 
  
  outcome <- c(a,b,c,d,e,f,g,h,i,j)
  stats <- data.frame(metric, outcome)
  
  # Export overall stats file
  write.table(stats, paste0(output, "_stats.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  message("Exported data")
  
  bam_primary <- bam_primary %>% 
    dplyr::mutate(above=coverage>0.95)
  
  # Histogram of coverage
  pdf(paste0(output, "_coverage_fraction.pdf"), width=6, height=6)
  plot1 <- ggplot(data=bam_primary, aes(x=coverage, fill=above)) +
    geom_histogram(bins = 180, show.legend = FALSE) +
    geom_vline(aes(xintercept=0.95), color="black", linetype="dashed", size=0.5) +
    xlim(0.5,1) +
    theme_classic(base_size=16) +
    xlab("Coverage Fraction") +
    ylab("Count") +
    scale_fill_manual(values = c("gray", "steelblue3"))
  print(plot1)
  dev.off()
  
  # Histogram of coverage vs length
  pdf(paste0(output, "_density.pdf"), width=8, height=5)
  plot2 <- ggplot() +
    geom_hex(data=bam_primary, aes(x=seqlengths, y=coverage, fill = stat(log(count))), bins=100) +
    stat_smooth(data=bam_primary, aes(x=seqlengths, y=coverage), color="lavender", se=TRUE, size=0.5, level=0.95) +
    xlim(0,15000) +
    ylim(0,1) +
    xlab("Known Transcript Length") +
    ylab("Coverage Fraction") +
    scale_fill_viridis_c() + 
    theme_classic(base_size=16)
  print(plot2)
  dev.off()
  
  # Secondary alignments bar chart
  pdf(paste0(output, "_sec_alns.pdf"), width=8, height=5)
  plot3 <- ggplot(alignments) +
    geom_bar(stat='identity', aes(x=nbrSecondary, y=prop), fill = "steelblue3") +
    xlab("Number of Secondary Alignments") +
    ylab("Proportion of Reads") +
    ylim(0,1) +
    scale_x_discrete(breaks = alignments$nbrSecondary, labels = alignments$nbrSecondary) +
    theme_classic(base_size=16) 
  print(plot3)
  dev.off() 
  
  # Histogram unqiue transcript lengths
  pdf(paste0(output, "_transcript_length_distribution.pdf"), width=6, height=6)
  plot4 <- ggplot(data=length_per_unique_transcript, aes(x=seqlengths)) +
    geom_histogram(bins = 180, show.legend = FALSE, fill="steelblue3") +
    theme_classic(base_size=16) +
    xlim(0,10000) +
    xlab("Known Transcript Length") +
    ylab("Count")
  print(plot4)
  dev.off()
  
  message("Complete")
  
}

suppressWarnings(
main())
