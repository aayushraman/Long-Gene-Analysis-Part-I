#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 19th April 2018
#
# Program is used for:
# 1. Analysis of NVS SEQC Microarray Dataset
#############################################################################

SEQC.microarray.analysis <- function(){
  #cat("Running SEQC Microarray Dataset\n")
  load(file = "../dat-infofiles/hg19.GPL17930.RData")
  
  ## Details about the experiment
  experiment <- pmid2MIAME("25150838")
  experiment@name <- "SEQC/MAQC-III Consortium"
  experiment@lab <- "SEQC/MAQC-III Consortium"
  experiment@url <- "https://www.ncbi.nlm.nih.gov/pubmed/25150838"
  annotation <- "HuGene2/GPL17930/Affymetrix Feature Matrix"
  dat.norm <- read.table(file = 
                        "../dat/GEO/GSE56457_SEQC_USF_GEOSub_PrimeView.txt", 
                         sep = "\t", stringsAsFactors = FALSE, quote = "", 
                         header = TRUE, na.strings = FALSE)
  A.samples <- dat.norm %>% dplyr::select(contains("_A_")) %>% colnames()
  B.samples <- dat.norm %>% dplyr::select(contains("_B_")) %>% colnames()
  C.samples <- dat.norm %>% dplyr::select(contains("_C_")) %>% colnames()
  D.samples <- dat.norm %>% dplyr::select(contains("_D_")) %>% colnames()
  genotypes <- factor(c(rep("A", 4), rep("B",4), rep("C",4),rep("D",4)),
                      levels = c("A", "B", "C", "D"))
  sample.info <- factor(c(rep("Human Ref.", 4), rep("Brain", 4),
                          rep("3Ref:1Brain", 4), rep("1Ref:3Brain",4)),
                levels = c("Human Ref.", "Brain", "3Ref:1Brain","1Ref:3Brain"))
  MDSplot.array <- MDSplot(data = 2^dat.norm[,c(2:17)], genotypes = genotypes)
  dat.norm <- aggregate(. ~ GeneSymbol, data = dat.norm, mean)
  rownames(dat.norm) <- dat.norm$GeneSymbol
  featureData <- read.table(file = 
            "../dat-infofiles/annot.files/RefSeq.hg19.Gene.TransciptLength.txt", 
                            sep = "\t", header = TRUE, quote = "", fill = TRUE,
                            stringsAsFactors = FALSE, na.strings = FALSE)
  sum(dat.norm$GeneSymbol %in% featureData$gene.name)
  hg19.gene.annot <- aggregate(. ~ gene.name, data = 
                    unique(featureData[,c("gene.name","gene.length")]), max)
  hg19.tx.annot <- aggregate(. ~ gene.name, data = 
                unique(featureData[,c("gene.name","transcript.length")]), max)
  dat.norm[,c(2:17)] <- data.frame(2^dat.norm[,c(2:17)])
  dat.annot <- inner_join(x = dat.norm, y = hg19.gene.annot, by = 
                              c("GeneSymbol" = "gene.name"))
  B1.samples <- B.samples[1:2]
  B2.samples <- B.samples[3:4]
  ratio.B1.B2.gl <- logofMeans.between.A.B(dat = dat.annot, B1.samples,
                                           B2.samples)
  Complot.B.array <- gabels.plot(mat = ratio.B1.B2.gl[,c("logFC.crude", 
                                                         "gene.length")]) 
  A1.samples <- A.samples[1:4]
  A2.samples <- A.samples[2:3]
  ratio.A1.A2.gl <- logofMeans.between.A.B(dat = dat.annot, A2.samples, 
                                           A1.samples)
  Complot.A.array <- gabels.plot(mat = ratio.A1.A2.gl[,c("logFC.crude", 
                                                         "gene.length")]) 
  ratio.ABC.gl <- logofMeans.between.ABC(dat = dat.annot, A.samples, B.samples, 
                                         C.samples)
  Complot.Beta.array <- gabels.plot(mat = ratio.ABC.gl[,c("logFC.crude",
                                                          "gene.length")],
                                    y.axis = "Mean Log2 (ß ratio)")
  Complot.Beta.array <- Complot.Beta.array + coord_cartesian(ylim = c(1.05,2.3))
  array.ABC.gl <- ratio.ABC.gl[,c("GeneSymbol", "logFC.crude", "gene.length")]
  save(file = "../dat/SEQC/SEQC_Microarray.RData", array.ABC.gl)
  ratio.ABC.tx <- logofMeans.between.ABC(dat = dat.annot, A.samples, B.samples,
                                         C.samples)
  ratio.ABC.tx <- inner_join(x = ratio.ABC.tx, y = hg19.tx.annot, 
                             by = c("GeneSymbol" = "gene.name"))
  array.ABC.tl <- gabels.plot(mat = unique(ratio.ABC.tx[,c("logFC.crude",
                                                         "transcript.length")]), 
                              length.type = "Transcript", 
                              y.axis = "Mean Log2 (ß ratio)")
  figures.array <- list(figure4.array = list(plotC = MDSplot.array, 
                                             plotD = Complot.Beta.array),
                        figureS6.array = list(plotB1 = Complot.B.array, 
                                              plotB2 = Complot.A.array,
                                              plotB3 = array.ABC.tl))
  return(figures.array)
}