# Replace sample1 with your file name or sample name
# Usage: Rscript BamSlam.R rna/cdna yourfile.bam ouputprefix

# GenomicAlignments/Features package is from bioconductor, need to install bioconductor then run:
# BiocManager::install("GenomicFeatures")

# Function for importing BAM file
import_bam_file <- function(bamfile, type) {
  bam <- GenomicAlignments::readGAlignments(bamfile, 
                                            use.names = TRUE,
                                            param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                                                 what = c("qname","flag","mapq")))
  
  message("Imported bam file")
  
  # Summarise CIGAR strings
  cigar_table <- cigarOpTable(bam@cigar)
  
  # CIGAR types 
  col_names_extract <- c("M", "I", "D", "S", "H")
  
  # Add summarised CIGAR strings
  mcols(bam)[paste0("nbr", col_names_extract)] <- mapply(function(col) cigar_table[, col], col_names_extract)
  
  message("Summarised CIGAR strings")
  
  bam_1 <- as.data.table(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::select(-cigar, -njunc)
  
  if (type == "cdna") {
    bam_2 <- bam_1[flag %in% c(0,16,256,272), ]
  } else if (type == "rna") {
    bam_2 <- bam_1[flag %in% c(0,256), ]
  } else {
    message("Sequencing type missing. Please enter either: cdna rna")
  }
  
  # Extract known mapped isoform lengths
  lengths <- as.data.frame(bam@seqinfo) %>% 
    dplyr::select(-isCircular, -genome) %>% 
    tibble::rownames_to_column("seqnames")
  
  # Merge known lengths and calculate transcript coverage
  bam_data <- left_join(bam_2, lengths, by="seqnames")

  message("Extracted known transcript lengths")
  
  return(bam_data)
  
}

# Function for calculating read coverages
get_read_coverages <- function(bam_data, output) {

  bam_data <- bam_data %>% 
    dplyr::mutate(aligned_length = nbrM + nbrI,
                  read_length = nbrS + nbrH + nbrM + nbrI,
                  aligned_fraction = aligned_length/read_length,
                  read_accuracy=(nbrM+nbrI+nbrD-NM)/(nbrM+nbrI+nbrD),
                  read_coverage=as.numeric(width/seqlengths)) %>% 
    tidyr::drop_na(read_coverage) %>% 
    rename(read_id = qname,
           transcript_id = seqnames,
           transcript_length = seqlengths) %>% 
    group_by(read_id) %>% 
    mutate(num_secondary_alns = n()-1) %>% ungroup()
  
  # get secondary alignments data
  alignments <- bam_data %>% 
    group_by(read_id) %>% 
    slice_head(n=1) %>% ungroup() %>% 
    group_by(num_secondary_alns) %>% 
    summarise(total = n()) %>%
    mutate(prop = total / sum(total))
  
  message("Calculated read coverages") 
  
  bam_data$num_secondary_alns <- as.factor(bam_data$num_secondary_alns)
  alignments$num_secondary_alns <- as.factor(alignments$num_secondary_alns)
  
  bam_export <- bam_data %>% dplyr::select(read_id, transcript_id, start, end, flag, mapq, AS, tp, aligned_length, 
                                           read_length, aligned_fraction, read_accuracy, transcript_length, read_coverage, num_secondary_alns)
  write.csv(bam_export, file = paste0(output, "_data.csv"), quote=F, row.names=F) 
  
  return(list(bam_data, alignments))
}

# Function for selecting the best AS per read
get_alignment_data_and_summarise <- function(bam_data, output) {
  bam_primary <- bam_data %>% 
    group_by(read_id) %>% 
    # arrange so P is before S
    arrange(tp) %>% 
    arrange(read_id, desc(AS)) %>% 
    dplyr::slice_head(n=1)
  
  bam_per_unique_transcript <- bam_primary %>% 
    dplyr::group_by(transcript_id) %>% 
    summarise(read_coverage = median(read_coverage, na.rm = TRUE),
              transcript_length = max(transcript_length))
  
  write.csv(bam_per_unique_transcript, file = paste0(output, "_transcript_level_data.csv"), quote=F, row.names=F)  
  
  # begin stats summary creation
  a <- nrow(bam_primary)
  b <- sum(bam_primary$read_coverage > 0.95)
  c <- b/nrow(bam_primary)*100
  d <- median(bam_primary$read_coverage)
  e <- median(bam_primary$aligned_length)
  f <- median(bam_primary$read_accuracy)*100
  g <- sum(bam_data$num_secondary_alns == 0)
  h <- g/nrow(bam_primary)*100
  i <- nrow(bam_per_unique_transcript)
  j <- median(bam_per_unique_transcript$read_coverage)
  k <- median(bam_per_unique_transcript$transcript_length)
  
  # Make stats
  metric <- c(
    "Sample",
    "Total number of reads",
    "Number of reads representing full-length transcripts",
    "Percentage of reads representing full-length transcripts",
    "Median coverage fraction of transcripts (primary alignments)",
    "Median alignment length (primary alignments)",
    "Median accuracy (primary alignments)",
    "Number of reads with no secondary alignments",
    "Percentage of reads with no secondary alignments",
    "Number of unique transcripts identified",
    "Median coverage fraction of all unique transcripts",
    "Median length of all unique transcripts identified") 
  
  outcome <- c(output,a,b,c,d,e,f,g,h,i,j,k)
  stats <- data.frame(metric, outcome)
  
  # Export overall stats file
  write.table(stats, file = paste0(output,"_stats.csv"), sep=",", quote=F, col.names = FALSE, row.names = FALSE) 
  message("Exported data")
  
  return(bam_primary)

}

# Function for plotting coverage
plot_coverage <- function(bam_primary, output) {
  bam_primary <- bam_primary %>% 
    dplyr::mutate(above=read_coverage>0.95)
  
  pdf(paste0(output, "_coverage_fraction.pdf"), width=6, height=6)
  plot <- ggplot(data=bam_primary, aes(x=read_coverage, fill=above)) +
    geom_histogram(bins = 180, show.legend = FALSE) +
    geom_vline(aes(xintercept=0.95), color="black", linetype="dashed", linewidth=0.5) +
    xlim(0.5,1) +
    theme_bw(base_size=14) +
    xlab("Coverage fraction") +
    ylab("Number of reads") +
    scale_fill_manual(values = c("gray", "steelblue3"))
  suppressMessages(print(plot))
  dev.off()
}

# Function for plotting coverage vs length
plot_coverage_vs_length <- function(bam_primary, output) {
  pdf(paste0(output, "_density.pdf"), width=8, height=5)
  plot <- ggplot() +
    geom_hex(data=bam_primary, aes(x=transcript_length, y=read_coverage, fill = after_stat(log(count))), bins=100) +
    stat_smooth(data=bam_primary, aes(x=transcript_length, y=read_coverage), color="lavender", se=TRUE, size=0.5, level=0.95) +
    xlim(0,15000) +
    ylim(0,1) +
    xlab("Known transcript length (nt)") +
    ylab("Coverage fraction") +
    scale_fill_viridis_c() + 
    theme_bw(base_size=14)
  suppressMessages(print(plot))
  dev.off()
}

# Function for plotting secondary alignments
plot_secondary_alignments <- function(alignments, output) {
  pdf(paste0(output, "_sec_alns.pdf"), width=8, height=5)
  plot <- ggplot(alignments) +
    geom_bar(stat='identity', aes(x=num_secondary_alns, y=prop), fill = "steelblue3") +
    xlab("Number of secondary alignments") +
    ylab("Proportion of reads") +
    ylim(0,1) +
    scale_x_discrete(breaks = alignments$num_secondary_alns, labels = alignments$num_secondary_alns) +
    theme_bw(base_size=14) 
  suppressMessages(print(plot))
  dev.off() 
}

# Function for plotting transcript length distribution
plot_transcript_length_distribution <- function(bam_primary, output) {
  
  length_per_unique_transcript <- bam_primary %>% 
    dplyr::group_by(transcript_id) %>% 
    slice_head(n=1)
  
  pdf(paste0(output, "_transcript_length_distribution_per_distinct_transcript.pdf"), width=6, height=6)
  plot <- ggplot(data=length_per_unique_transcript, aes(x=transcript_length)) +
    geom_histogram(bins = 180, show.legend = FALSE, fill="steelblue3") +
    theme_bw(base_size=14) +
    xlim(0,10000) +
    xlab("Known transcript length (nt)") +
    ylab("Transcript count")
  suppressMessages(print(plot))
  dev.off()
  
  pdf(paste0(output, "_transcript_length_distribution_per_read.pdf"), width=6, height=6)
  plot <- ggplot(data=bam_primary, aes(x=transcript_length)) +
    geom_histogram(bins = 180, show.legend = FALSE, fill="steelblue3") +
    theme_bw(base_size=14) +
    xlim(0,10000) +
    xlab("Known transcript length (nt)") +
    ylab("Read count")
  suppressMessages(print(plot))
  dev.off()
  
}

# Function for plotting accuracy
plot_accuracy <- function(bam_primary, output) {
  pdf(paste0(output, "_read_accuracy.pdf"), width=6, height=6)
  plot <- ggplot(data=bam_primary, aes(x=read_accuracy, y=after_stat(scaled))) +
    geom_density(alpha = 0.4, show.legend = FALSE, fill="steelblue3") +
    theme_bw(base_size=14) +
    xlim(0.5,1) +
    xlab("Read accuracy") +
    ylab("Density")
  suppressMessages(print(plot))
  dev.off()
}

# main function
main <- function() {
 
  # inputs from command line
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
    library(hexbin)
    library(data.table)
  })
  
  options(dplyr.summarise.inform = FALSE)
  
  # call functions
  bam_tidy <- import_bam_file(bamfile, type)
  coverages_list_data <- get_read_coverages(bam_tidy, output)
  
  bam_data <- coverages_list_data[[1]]
  alignments <- coverages_list_data[[2]]
  
  bam_primary <- get_alignment_data_and_summarise(bam_data, output)
  
  # plotting functions
  plot_coverage(bam_data, output)
  plot_coverage_vs_length(bam_data, output)
  plot_secondary_alignments(alignments, output)
  plot_transcript_length_distribution(bam_primary, output)
  plot_accuracy(bam_data, output)
  
  message("Complete")
  
}

suppressWarnings(
  main())
