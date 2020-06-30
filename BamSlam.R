
# Replace sample1 with your file name or sample name
# Usage: Rscript BamSlam.R yourfile.bam gencode.gtf ouputprefix

# GenomicAlignments/Features package is from bioconductor, need to install bioconductor then run:
# BiocManager::install("GenomicFeatures")

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  bamfile <- args[1]
  gtffile <- args[2]
  output <- args[3]
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(GenomicFeatures)
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(tidyr)
    library(MASS)
    library(viridis)
  })
  
  # Import bam
  bam <- readGAlignments(bamfile, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                              what = c("qname","flag", "rname", 
                                                       "pos", "mapq", "seq", "qual")))
  # Expand cigar strings
  ops <- GenomicAlignments::CIGAR_OPS
  wdths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = ops)
  keep.ops <- GenomicAlignments::explodeCigarOps(cigar(bam), ops = ops)
  explodedcigars <- IRanges::CharacterList(relist(paste0(unlist(wdths), 
                                                         unlist(keep.ops)), wdths))
  for (opts in setdiff(GenomicAlignments::CIGAR_OPS, "=")) {
    mcols(bam)[[paste0("nbr", opts)]] <- 
      sapply(explodedcigars, function(cg) sum(as.numeric(gsub(paste0(opts, "$"), "", cg)), na.rm = TRUE))
  }
  
  # Add columns of interest
  bam_data <- data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::mutate(alignedLength = nbrM + nbrI) %>% 
    dplyr::mutate(readLength = nbrS + nbrH + nbrM + nbrI) %>% 
    dplyr::mutate(alignedFraction = alignedLength/readLength) %>% 
    dplyr::mutate(accuracy=(nbrM+nbrI+nbrD-NM)/(nbrM+nbrI+nbrD))
  
  # Make the txdb
  txs <- makeTxDbFromGFF(gtffile, format="gtf")
  txLengths <- transcriptLengths(txs, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  lengths <- data.frame(txLengths$tx_name, txLengths$tx_len)
  bam_data$seqnames <- as.character(bam_data$seqnames)
  
  # Make it work with alignments to gencode fasta
  if((grepl("\\|", bam_data$seqnames)) == TRUE) {
    bam_data <- data.frame(bam_data, do.call(rbind, strsplit(bam_data$seqnames, "\\|"))[,1])
    bam_data$do.call.rbind..strsplit.bam_data.seqnames............1. <- as.character(bam_data$do.call.rbind..strsplit.bam_data.seqnames............1.)
    bam_data$transcript <- bam_data$do.call.rbind..strsplit.bam_data.seqnames............1.
    bam_data$seqnames <- NULL
    bam_data$do.call.rbind..strsplit.bam_data.seqnames............1. <- NULL
  } else {
    bam_data$transcript <- bam_data$seqnames
    bam_data$seqnames <- NULL
  }
  
  # Merge known tx length
  bam_data <- merge(bam_data, lengths, by.x="transcript", by.y="txLengths.tx_name", all.x=TRUE)
  
  bam_data <- bam_data %>% 
    dplyr::mutate(coverage=width/txLengths.tx_len)
  
  # Filter 
  bam_data <- bam_data %>% 
    dplyr::filter(strand == "+")
  
  #bam_filtered <- bam_data %>% 
    #dplyr::filter(end > (txLengths.tx_len - 100))
  
  bam_filtered <- bam_data %>% 
    group_by(qname) %>% 
    arrange(qname, desc(AS)) %>% 
    filter(AS / max(AS) >= 0.9)
  
  bam_filtered <- bam_filtered %>% 
    group_by(qname) %>% 
    arrange(qname, desc(alignedFraction)) %>% 
    filter(alignedFraction / max(alignedFraction) >= 0.9)
  
  bam_filtered <- bam_filtered %>% 
    group_by(qname) %>% 
    arrange(qname, desc(accuracy)) %>% 
    filter(accuracy / max(accuracy) >= 0.9)

  #issue
  # Export files
  write.table(bam_data, file = paste0(output, "_data.txt"), sep="\t", quote=F, col.names = T, row.names=F) 

  bam_primary <- subset(bam_data, tp == "P")
  a <- sum(bam_primary$coverage > 0.95)
  b <- nrow(bam_primary)
  c <- a/nrow(bam_primary)*100
  d <- median(bam_primary$coverage)
  e <- median(bam_primary$accuracy)*100
  f <- length(unique(bam_primary$qname))
  g <- f/nrow(bam_primary)*100
  
  # Make stats
  metric <- c(
    "Number of reads representing full-length transcripts:",
    "Out of total primary alignments:",
    "Percentage of reads representing full-length transcripts:",
    "Median coverage fraction of transcripts:",
    "Median accuracy of primary alignments:",
    "Number of reads that identify unique transcripts:",
    "Percentage of reads that identify unique transcripts:") 
  
  outcome <- c(a,b,c,d,e,f,g)
  stats <- data.frame(metric, outcome)
  
  # Export overall stats file
  write.table(stats, paste0(output, "_stats.txt"), sep="\t", quote=F, row.names=F, col.names=F)
 
  bam_primary <- bam_primary %>% 
    dplyr::mutate(above=coverage>0.95)
  
  # Histogram
  pdf(paste0(output, "_coverage_fraction.pdf"), width=6, height=6)
  ggplot(data=bam_primary, aes(x=coverage, fill=above)) +
    geom_histogram(bins = 180, show.legend = FALSE) +
    geom_vline(aes(xintercept=0.95), color="black", linetype="dashed", size=0.5) +
    xlim(0.5,1) +
    theme_classic(base_size=16) +
    xlab("Coverage Fraction") +
    ylab("Count") +
    scale_fill_manual(values = c("gray", "steelblue3"))
  dev.off()
  
  # 2D Density
  pdf(paste0(output, "_2d_density.pdf"), width=8, height=5)
  ggplot() + 
    stat_density_2d(data=bam_primary, aes(x=txLengths.tx_len, y=coverage, fill = ..level..), geom = "polygon",
                    position = "identity",
                    na.rm = TRUE,
                    show.legend = TRUE,
                    inherit.aes = TRUE,
                    bins = 50
    ) +
    stat_smooth(data=bam_primary, aes(x=txLengths.tx_len, y=coverage), color="lavender", se=TRUE, size=1, level=0.95) +
    xlim(0,8000) +
    ylim(0,1) +
    xlab("Known Transcript Length") +
    ylab("Coverage Fraction") +
    scale_fill_viridis_c() +
    theme_classic(base_size=16)
  dev.off() 
  
}

main()

