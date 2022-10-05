#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 14th April 2018
#
# Program is used for:
# 1. Long and Short Gene analysis for Gabel's VC Dataset
#############################################################################

## dataset
gabel.analysis <- function(){
  gabel.dataset <- read.table("../dat/Gabel_VC_Dataset/GBCounts.txt",sep = "\t", 
                              stringsAsFactors=F, header = TRUE, row.names = 1)
  head(gabel.dataset)
  gabel.DESeq <- read.table("../dat/Gabel_VC_Dataset/GB_DESeq2.txt",sep = "\t", 
                            stringsAsFactors=F, header = TRUE, row.names = 1)
  genotype <- data.frame(c(rep("WT", 3), rep("KO", 3)))
  rownames(genotype) = colnames(gabel.dataset)
  colnames(genotype) = "genotypes"
  genotype$genotypes <- relevel(genotype$genotypes, "WT")
  dds <- DESeqDataSetFromMatrix(gabel.dataset, genotype, design = ~ genotypes)
  dds <- estimateSizeFactors(dds)
  dat <- counts(dds, normalized = TRUE)
  ## rlog transformation
  rld <- rlogTransformation(dds, blind = TRUE)
  pca.rld = plotPCA(rld, intgroup = c("genotypes"), returnData = TRUE)
  percentVar = round(attr(pca.rld, "percentVar")*100 ,2)
  title <- "PCA Plot for KO Dataset"
  #print(plotPCA(rld, intgroup = c("genotypes")))
  #MDSplot(data = assay(rld), genotypes = genotype)
  
  ## MECP2 conc. and box plot
  ind <- which(rownames(dat) == "Mecp2")  ## Mecp2, "NM_001110792"
  plot.dat <- data.frame(labels = factor(colnames(dat[,c(4:6,1:3)]), 
                                         levels = colnames(dat[,c(4:6,1:3)])),
                        Normalized.Counts = as.vector(
                            as.matrix(log2(dat[ind,c(4:6,1:3)])+1)), 
                        genotypes = relevel(genotype$genotypes[c(4:6,1:3)],
                                            "KO"))
  ## merge the DESeq and counts file 
  gabel.DESeq <- merge(x = gabel.dataset, y = gabel.DESeq, by = "row.names")
  rownames(gabel.DESeq) <- gabel.DESeq$Row.names
  gabel.DESeq <- gabel.DESeq[,-1]
  #head(gabel.DESeq)
  
  ## Expression Test
  dds <- DESeq(dds)
  res <- results(dds)
  res$norm.counts <- counts(dds, normalized=TRUE)
  #sum(res$padj < 0.05, na.rm = TRUE)
  
  ## other data frame for the plots
  dat2 <- as.data.frame(res)
  dat2$logFC.crude <- apply(dat2[,c(7:12)], 1, 
                        function(r){log2((mean(r[4:6])+1)/(mean(r[1:3])+1))})
  
  ## reading the annotation file
  mm10.annot <- read.table(
      file = "../dat-infofiles/annot.files/mm10_ucsc.annot.txt", header = FALSE, 
      sep = "\t", stringsAsFactors = FALSE, quote = "", na.strings = FALSE)
  mm10.annot <- mm10.annot[,c(1,3:6)]
  colnames(mm10.annot) <- c("Gene_name", "chr", "strand", "tx_start", "tx_end")
  mm10.annot$gene.length <- mm10.annot$tx_end - mm10.annot$tx_start + 1
  mm10.annot <- aggregate(.~ Gene_name, data = mm10.annot[,c("Gene_name", 
                                                             "gene.length")],max)
  
  ## merge degs dat with mm10
  dat3 <- merge(dat2, mm10.annot, by.x = "row.names", by.y = "Gene_name")
  #p1 <- gabels.plot(mat = dat3[,c("logFC.crude", "gene.length")], comp.between = "(KO/WT)")
  #p1$plot + coord_cartesian(ylim = c(-0.15,0.15))
  gene.type <- ifelse(dat3$padj < 0.05 & !is.na(dat3$padj), 
                     ifelse(dat3$gene.length > 100e3, "Long Genes", 
                            "Short Genes"),"Not Stat. Signif.")
  # Plot.Scatter(dat = dat3[,c("Row.names", "logFC.crude", "padj", "gene.length")], 
  #              log2FC = 0, comp.between = "(KO/WT)", pval = 0.05)
  ## gabel plot for WT/WT and KO/WT for Visual Cortex
  colnames(dat3)[1] <- "gene.name"
  dat3$logFC.WT <- apply(X = dat3[,c(8:10)], 1, function(r) (log2((r[3]+1)/
                                                                (r[c(1)]+1))))
  p1 <- overlay.gabels.plot(mat = dat3[,c(16,14,15)], comp.between1 = "(WT/WT)", 
                            comp.between2 = "(KO/WT)")
  cat("\n\n Printing Supp. Figure 2 -- Visual Cortex (RNA-seq)\n\n")
  print(plot_grid(p1$plot1 + coord_cartesian(ylim = c(-0.2,0.2)), p1$plot2 + 
                      coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
}