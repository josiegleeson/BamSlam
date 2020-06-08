# Written by Josie Gleeson 2020 
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
  })
  
  # Import bam
  bam <- readGAlignments(bamfile, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                              what = c("qname","flag", "rname", 
                                                       "pos", "mapq", "seq", "qual")))
  # Start building csv data
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
  txs <- makeTxDbFromGFF("gencode.v31.annotation.gtf", format="gtf")
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
  bam_filtered <- bam_data %>% 
    dplyr::filter(strand == "+")
  
  bam_filtered <- bam_filtered %>% 
    group_by(qname) %>% 
    dplyr::filter((end <= txLengths.tx_len) &&  (end > txLengths.tx_len - 100)) 
  
  bam_filtered <- bam_filtered %>% 
    group_by(qname) %>% 
    arrange(qname, desc(AS)) %>% 
    filter(AS / max(AS) >= 0.98)
  
  bam_filtered <- bam_filtered %>% 
    group_by(qname) %>% 
    arrange(qname, desc(alignedFraction)) %>% 
    filter(alignedFraction / max(alignedFraction) >= 0.9)
  
  bam_filtered <- bam_filtered %>% 
    group_by(qname) %>% 
    arrange(qname, desc(accuracy)) %>% 
    filter(accuracy / max(accuracy) >= 0.9)
  
  # Export file as a sam file
  bam_export <- subset(bam_filtered, select=c("qname", "flag", "transcript", "start", "mapq", "cigar", "seq", "qual"))
  bam_export <- bam_export %>% 
    dplyr::mutate(rnext="*", pnext=0, tlen=0)
  
  col_order <- c("qname", "flag", "transcript", "start", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual")
  bam_export <- bam_export[, col_order]
  
  bam_export$seq[bam_export$seq==""]<-"*"
  bam_export$qual[bam_export$qual==""]<-"*"
  
  # Export files
  write.csv(bam_data, file = paste0(output, "_data.csv"), quote=F, col.names = F) 
  write.table(bam_export, file = paste0(output, "_pre_filtered.sam"), sep="\t", quote=F, col.names = F, row.names = F)

  bam_primary <- subset(bam_data, tp == "P")
  a <- sum(bam_primary$coverage > 0.95)
  b <- nrow(bam_primary)
  c <- a/nrow(bam_primary)*100
  d <- median(bam_primary$coverage)
  e <- median(bam_primary$accuracy)*100
  f <- nrow(bam_data)
  g <- nrow(subset(bam_data, tp == "P"))
  h <- nrow(subset(bam_data, tp == "S"))
  i <- nrow(bam_filtered)
  j <- nrow(subset(bam_filtered, tp == "P"))
  k <- nrow(subset(bam_filtered, tp == "S"))
  
  # Make stats
  metric <- c(
    "Number of reads representing full-length transcripts:",
    "Out of total primary alignments before filtering:",
    "Percentage of reads representing full-length transcripts:",
    "Median coverage fraction of transcripts:",
    "Median accuracy of primary alignments:",
    "Alignments before filtering:",
    "Number of initial alignments that are primary:",
    "Number of initial alignments that are secondary:",
    "Alignments left after filtering:",
    "Number of filtered alignments that are primary:",
    "Number of filtered alignments that are secondary:")
  outcome <- c(a,b,c,d,e,f,g,h,i,j,k)
  stats <- data.frame(metric, outcome)
  
  # Export overall stats file
  write.table(stats, paste0(output, "_stats.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  
}

main()
