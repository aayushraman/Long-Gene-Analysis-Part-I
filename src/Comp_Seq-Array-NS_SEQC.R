#########################################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 18th April 2018
#
# Program is used for:
# 1. RTT nCounter data analysis -- Comparison of log2FC.crude vs DESeq2 log2FoldChange
#########################################################################################

Complog2FC_tech <- function(){
  
    ## SEQC RNA-Seq and Microarray
    load("../dat/SEQC/SEQC_RNA-Seq.RData")
    load("../dat/SEQC/SEQC_Microarray.RData")
  
    ## reading human nanoString Dataset
    seqc.NanoString <- read.table("../dat/SEQC/nCounterDataset.txt", sep = "\t", 
                                  stringsAsFactors = FALSE, quote = "", 
                                  header = TRUE, na.strings = FALSE)
    seqc.NanoString.norm <- NanoStringNorm(x = seqc.NanoString, anno = NA,
                                      CodeCount = 'sum', Background = 'mean',
                                      SampleContent = 'housekeeping.sum',
                                      round.values = FALSE,take.log = FALSE,
                                      return.matrix.of.endogenous.probes = TRUE)
    genotypes <- factor(rep(c(rep("A",3),rep("B",3),rep("C",3),rep("D",3)),2), 
                      levels = c("A", "B", "C","D"))
    #boxPlot(data = log2(seqc.NanoString.norm+1), title = 
    # "Log Transformed Nano string Dataset", samples = genotypes)
    p1 <- MDSplot(data = log2(seqc.NanoString.norm+1), genotypes = genotypes)
    seqc.NanoString.norm <- data.frame(seqc.NanoString.norm)
    UHR <- seqc.NanoString.norm %>% dplyr::select(contains("UHR_")) %>% 
                                                                    colnames()
    HBR <- seqc.NanoString.norm %>% dplyr::select(contains("HBR_")) %>% 
                                                                    colnames()
    C.type <- seqc.NanoString.norm %>% dplyr::select(contains("C_")) %>% 
                                                                    colnames()
    seqc.NanoString.norm <- seqc.NanoString.norm[, order(genotypes)]
    genotypes <- genotypes[order(genotypes)]
    seqc.NanoString.norm <- as.data.frame(seqc.NanoString.norm)
    seqc.NanoString.norm$Mean.B <- apply(seqc.NanoString.norm[,HBR], 1, 
                                         function(r) {(mean(r))})
    seqc.NanoString.norm <- seqc.NanoString.norm %>% rownames_to_column()
    hg19.annot <- read.table(
    file = "../dat-infofiles/annot.files/RefSeq.hg19.Gene.TransciptLength.txt", 
    header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    hg19.annot <- hg19.annot[,c("gene.name", "gene.length")]
    hg19.annot <- aggregate( .~ gene.name, data = hg19.annot, max)
    seqc.NanoString.norm <- inner_join(x = seqc.NanoString.norm, y = hg19.annot, 
                                       by = c("rowname" = "gene.name"))
    colnames(seqc.NanoString.norm)[1] <- "gene.name"
    ind <- seqc.NanoString.norm$gene.length > 100e3
    cols <- ifelse(ind == TRUE, "blue", "red")
    sum(seqc.NanoString.norm$Mean.B > 256 & seqc.NanoString.norm$gene.length > 
                    100e3) ## 132 out 184 are expressed
    p2 <- qplot(x = gene.length/1000, y = log2(Mean.B+1), 
                data = seqc.NanoString.norm, col = cols) + 
              xlab("Gene Length in KB") + ylab("Mean Gene Expr. (Brain Samples)") + 
              scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,
                                                                   1000,10000))+
              theme(plot.title = element_text(size = 14, face = "bold"),
                    legend.position = "none", 
                    axis.title.y = element_text(size = 18, colour = "black", 
                                                face = "bold"),
                    axis.title.x = element_text(size = 18, colour = "black", 
                                                face = "bold"),
                    axis.text.y = element_text(size = 18, colour = "black", 
                                               face = "bold"),
                    axis.text.x = element_text(size = 18, colour = "black", 
                                               face = "bold"), 
                    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    nCounter.B1.B2.gl <- logofMeans.between.A.B(dat = seqc.NanoString.norm %>% 
                                                dplyr::select(matches("HBR|gene.")),
                                                A.samples = HBR[1:3], 
                                                B.samples = HBR[4:6])
    p3 <- scatter.lm(unique(nCounter.B1.B2.gl[,c("gene.name", "logFC.crude", 
                                             "gene.length")]))
    nCounter.ABC.gl <- logofMeans.between.ABC(dat = seqc.NanoString.norm %>%
                                            dplyr::select(matches("UHR|HBR|C_|gene.")), 
                                            A.samples = UHR, B.samples = HBR,
                                            C.samples = C.type)
    p4 <- scatter.lm(unique(nCounter.ABC.gl[,c("gene.name", "logFC.crude",
                                               "gene.length")]))
    seq.ABC.gl <- seq.ABC.gl[seq.ABC.gl$gene.name %in% nCounter.ABC.gl$gene.name,]
    array.ABC.gl <- array.ABC.gl[array.ABC.gl$GeneSymbol %in% 
                                       nCounter.ABC.gl$gene.name,]
    seq.ABC.gl$longGene <- ifelse(seq.ABC.gl$gene.length > 100e3, "LongGene",
                                  "ShortGene")
    seq.ABC.gl$longGene <- factor(seq.ABC.gl$longGene, levels = c("ShortGene",
                                                                  "LongGene"))
    p5 <- boxPlot.comp(dat = seq.ABC.gl, type = "RNA-Seq")
    array.ABC.gl$longGene <- ifelse(array.ABC.gl$gene.length > 100e3, 
                                    "LongGene","ShortGene")
    array.ABC.gl$longGene <- factor(array.ABC.gl$longGene, levels = 
                                        c("ShortGene", "LongGene"))
    p6 <- boxPlot.comp(dat = array.ABC.gl, type = "Array")
    nCounter.ABC.gl$longGene <- ifelse(nCounter.ABC.gl$gene.length > 100e3, 
                                       "LongGene", "ShortGene")
    nCounter.ABC.gl$longGene <- factor(nCounter.ABC.gl$longGene, 
                                       levels = c("ShortGene", "LongGene"))
    p7 <- boxPlot.comp(dat = nCounter.ABC.gl, type = "nanoString")
    figures.log2FCcomp <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, 
                               p6 = p6, p7 = p7)
    return(figures.log2FCcomp)
}