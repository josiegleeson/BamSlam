# Replace sample1 with your file name or sample name
# Usage: Rscript BamSlam.R yourfile.bam gencode.gtf ouputprefix

# GenomicAlignments/Features package is from bioconductor, need to install bioconductor then run:
# BiocManager::install("GenomicFeatures")

# Lines 28-66 are from parts of: https://github.com/csoneson/NativeRNAseqComplexTranscriptome (Copyright (c) 2019 Charlotte Soneson)


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
  message("Imported bam file")
  
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
  
  # Get number of secondary and supplementary alignments
  tmp <- data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::rename(nbrJunctions = njunc) %>%
    dplyr::select(-cigar) 
  
  tmp2 <- as.data.frame(table(names(subset(bam, flag %in% c(0, 16)))))
  if (nrow(tmp2) == 0) tmp2 <- data.frame(Var1 = tmp$qname[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp2 %>% dplyr::rename(qname = Var1, nbrPrimaryAlignments = Freq))
  
  tmp3 <- as.data.frame(table(names(subset(bam, flag %in% c(256, 272)))))
  if (nrow(tmp3) == 0) tmp3 <- data.frame(Var1 = tmp$qname[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp3 %>% dplyr::rename(qname = Var1, nbrSecondaryAlignments = Freq))
  
  tmp4 <- as.data.frame(table(names(subset(bam, flag %in% c(2048, 2064)))))
  if (nrow(tmp4) == 0) tmp4 <- data.frame(Var1 = tmp$qname[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp4 %>% dplyr::rename(qname = Var1, nbrSupplementaryAlignments = Freq))
  
  tmp <- tmp %>% dplyr::mutate(nbrSecondaryAlignments = replace(nbrSecondaryAlignments, 
                                                                is.na(nbrSecondaryAlignments), 0),
                               nbrSupplementaryAlignments = replace(nbrSupplementaryAlignments, 
                                                                    is.na(nbrSupplementaryAlignments), 0))
  message("Created CIGAR string columns")
  
  bam_data <- tmp %>%
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
  
  message("Imported GTF")           
  
  # Merge known tx length
  bam_data <- merge(bam_data, lengths, by.x="transcript", by.y="txLengths.tx_name", all.x=TRUE)
  
  bam_data <- bam_data %>% 
    dplyr::mutate(coverage=width/txLengths.tx_len)
  
  bam_data$coverage <- as.numeric(bam_data$coverage)
  bam_data <- bam_data %>% 
    tidyr::drop_na(coverage)           
  
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
  
  unique_reads <-  bam_filtered %>% 
    group_by(qname) %>% 
    dplyr::filter(n()==1)
  
  message("Calculated transcript coverages")            
  
  # Export files
  bam_export <- subset(bam_data, select=c("transcript", "qwidth", "start", "end", "width", "qname", "flag", "mapq", "NM", "AS", "tp", "nbrM", "nbrI", "nbrD", "nbrN", "nbrS", "alignedLength", "readLength", "alignedFraction", "accuracy", "txLengths.tx_len", "coverage", "nbrSecondaryAlignments", "nbrSupplementaryAlignments"))
  write.csv(bam_export, file = paste0(output, "_data.csv"), sep=",", quote=F, col.names = T, row.names=F) 
  message("Exported data as csv")
  bam_sec <- bam_data %>% 
    dplyr::group_by(qname) %>% 
    dplyr::arrange(qname, desc(nbrSecondaryAlignments)) %>% 
    dplyr::slice(n=1)
  
  unique_reads <-  filter(bam_sec, nbrSecondaryAlignments == 0)
  
  bam_sec <- bam_sec %>% 
    group_by(nbrSecondaryAlignments) %>% 
    summarise(total = n()) %>%
    mutate(prop = total / sum(total))
  
  bam_sec$nbrSecondaryAlignments <- as.factor(bam_sec$nbrSecondaryAlignments)
  
  bam_primary <- subset(bam_data, flag == 0)
  
  bam_per_unique_transcript <- bam_primary %>% 
    dplyr::group_by(transcript) %>% 
    summarise(coverage = median(coverage, na.rm = TRUE))
  
  write.csv(bam_per_unique_transcript, file = paste0(output, "_transcript_level_data.csv"), sep=",", quote=F, col.names = T, row.names=F)  
  
  a <- sum(bam_primary$coverage > 0.95)
  b <- nrow(bam_primary)
  c <- a/nrow(bam_primary)*100
  d <- median(bam_primary$coverage)
  e <- median(bam_primary$accuracy)*100
  f <- nrow(unique_reads)
  g <- f/nrow(bam_primary)*100
  h <- nrow(bam_per_unique_transcript)
  i <- median(bam_per_unique_transcript$coverage)
  
  # Make stats
  metric <- c(
    "Number of reads representing full-length transcripts:",
    "Out of total primary alignments:",
    "Percentage of reads representing full-length transcripts:",
    "Median coverage fraction of transcripts:",
    "Median accuracy of primary alignments:",
    "Number of reads with no secondary alignments:",
    "Percentage of reads with no secondary alignments:",
    "Number of unique transcripts identified:",
    "Median coverage per transcript:") 
  
  outcome <- c(a,b,c,d,e,f,g,h,i)
  stats <- data.frame(metric, outcome)
  
  # Export overall stats file
  write.table(stats, paste0(output, "_stats.txt"), sep="\t", quote=F, row.names=F, col.names=F)
  message("Exported stats file")
  bam_primary <- bam_primary %>% 
    dplyr::mutate(above=coverage>0.95)
  message("Creating plots")
  # Histogram
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
  
  # 2D Density
  pdf(paste0(output, "_2d_density.pdf"), width=8, height=5)
  plot2 <- ggplot() + 
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
  print(plot2)
  dev.off() 
  
  pdf(paste0(output, "_sec_alns.pdf"), width=8, height=5)
  plot3 <- ggplot(bam_sec) +
    geom_bar(stat='identity', aes(x=nbrSecondaryAlignments, y=prop), fill = "steelblue3") +
    xlab("Number of Secondary Alignments") +
    ylab("Proportion of Reads") +
    ylim(0,1) +
    theme_classic(base_size=16)
  print(plot3)
  dev.off() 
  message("Complete")
}

main()
