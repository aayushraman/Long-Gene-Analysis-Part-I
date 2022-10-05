## DESeq2 Function
DESeqCalculation <- function(seqcData, genotype){
  dds <- DESeqDataSetFromMatrix(seqcData, genotype, design = ~ genotypes)
  dds <- estimateSizeFactors(dds)
  dat <- counts(dds, normalized=TRUE)
  #message("Total Number of reads generated at Novartis for A, B and C samples only: ",sum(round(colSums(assay(dds)), 1)))
  
  ## MECP2 conc.
  #message("\n\n MeCP2 concentration\n\n")
  ind <- which(rownames(dat) == "NM_004992")  ## "NM_001110792"
  plot.dat <- data.frame(labels = factor(colnames(dat), levels = colnames(dat)),
                         Normalized.Counts = as.vector(as.matrix(dat[ind,])),
                         genotypes = genotype$genotypes)
  # print(ggplot(plot.dat, aes(x = labels, y = Normalized.Counts, fill = genotypes)) + 
  #         geom_bar(stat="identity", width=0.75, position = position_dodge(width=0.5)) +
  #         xlab("") + ylab("Normalized Counts") + ggtitle(paste("Barplot for Gene MeCP2")) + theme_bw() +
  #         theme(plot.title = element_text(size = 14, face = "bold"),
  #               legend.text = element_text(size = 10, face = "bold"),
  #               legend.title = element_text(colour = "black", size = 12),
  #               axis.title.y= element_text(size = 10, colour = "black", face = "bold"),
  #               axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
  #               axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  return(dat)
}