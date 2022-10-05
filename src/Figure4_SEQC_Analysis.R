#############################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 19th April 2018
#
# Program is used for:
# 1. Analysis of NVS SEQC RNA-Seq Dataset
#############################################################################

figure4 <- function(){
    # figures.seq <- figures_SEQC()
    
    ## Figure 4
    load(file = "../results/SEQC_Figures.RData")
    cat("\n\n Fig. 4(A) -- SEQC RNA-seq\n\n")
    print(figures.seq$figure4$plotA + geom_point(size = 5))
    print(figures.seq$figure4$plotB)
    
    cat("\n\n Fig. 4(B) -- SEQC Array\n\n")
    print(figures.seq$figure4$plotC + geom_point(size = 5))
    print(figures.seq$figure4$plotD)
    
    cat("\n\n Fig. 4(E) -- RNA-seq\n\n")
    print(figures.seq$figure4$plotE)
    
    cat("\n\n Fig. 4(F) -- Microarray\n\n")
    print(figures.seq$figure4$plotF)
    
    cat("\n\n Fig. 4(F) -- Nanostring\n\n")
    print(figures.seq$figure4$plotG)
}

figureS6 <- function(){
    load(file = "../results/SEQC_Figures.RData")
    cat("\n\n Supplementary Fig. 6(A) -- SEQC RNA-seq\n\n")
    print(figures.seq$figureS6.seq$plotA1)
    print(figures.seq$figureS6.seq$plotA2)
    print(figures.seq$figureS6.seq$plotA3)
    cat("\n\n Supplementary Fig. 6(B) -- SEQC Array\n\n")
    print(figures.seq$figureS6.seq$plotB1)
    print(figures.seq$figureS6.seq$plotB2)
    print(figures.seq$figureS6.seq$plotB3)
}

figureS7 <- function(){
    load(file = "../results/SEQC_Figures.RData")
    cat("\n\n Supplementary Fig. 7(A) -- Total Count\n\n")
    print(figures.seq$figureS7$plotA1)
    print(figures.seq$figureS7$plotA2)
    cat("\n\n Supplementary Fig. 7(B) -- TMM (edgeR)\n\n")
    print(figures.seq$figureS7$plotB1)
    print(figures.seq$figureS7$plotB2)
}

figureS8 <- function(){
    load(file = "../results/SEQC_Figures.RData")
    print(figures.seq$figureS8$plotA + geom_point(size = 5))
    print(figures.seq$figureS8$plotB)
    print(figures.seq$figureS8$plotC)
    print(figures.seq$figureS8$plotD)
}

figureS9 <- function(){
    load(file = "../results/SEQC_Figures.RData")
    cat("\n\n Supplementary Fig. 9(A) -- Technical Brain Replicates \n\n")
    print(figures.seq$figureS9$plotA)
    cat("\n\n Supplementary Fig. 9(B) -- Comparison of Lib1/Lib2 \n\n")
    print(figures.seq$figureS9$plotB)
    cat("\n\n Supplementary Fig. 9(C) -- Comparison of Lib3/Lib1 \n\n")
    print(figures.seq$figureS9$plotC)
    cat("\n\n Supplementary Fig. 9(D) -- DEGs Comparison of Lib1/Lib2 \n\n")
    print(figures.seq$figureS9$plotD1)
    print(figures.seq$figureS9$plotD2)
    cat("\n\n Supplementary Fig. 9(D) -- DEGs Comparison of Lib3/Lib2 \n\n")
    print(figures.seq$figureS9$plotD3)
    print(figures.seq$figureS9$plotD4)
}

figures_SEQC <- function(){
    
    ## NVS dataset and hg19 annot
    counts.table <- read.table(file = "../dat/SEQC/GSE47774_SEQC_ILM_NVS.txt", 
                               header = TRUE, sep = "\t", quote = "", 
                               row.names = 1, stringsAsFactors = FALSE)
    hg19.annot <- read.table(file = 
            "../dat-infofiles/annot.files/RefSeq.hg19.Gene.TransciptLength.txt",
            header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    erccindex <- grep("ERCC",rownames(counts.table))
    counts.table <- counts.table[-erccindex,]
  
    ## sample types
    samples.type <- gsub(pattern = 
        "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_(.*)_(.*)", 
                         replacement = "\\3", x = colnames(counts.table))
    samples.barcode <- gsub(pattern = 
        "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_([ATGC]+)_(.*)", 
                            replacement = "\\6", x = colnames(counts.table))
    samples.libID <- gsub(pattern = 
        "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_([ATGC]+)_(.*)", 
                          replacement = "\\4", x = colnames(counts.table))
    samples.flowcell <- gsub(pattern = 
        "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_([ATGC]+)_(.*)",
                             replacement = "\\7", x = colnames(counts.table))
  
    ## Arguments for DESeq2 and its run
    ind <- which(samples.type == "A" | samples.type == "B" | 
                     samples.type == "C" | samples.type == "D")
    counts.table <- counts.table[,ind]
    colData <- as.data.frame(samples.type[ind])
    rownames(colData) <- colnames(counts.table)
    colnames(colData) <- "genotypes"
    exprs.dat <- DESeqCalculation(seqcData = counts.table, genotype = colData)
    MDSplot.seq <- MDSplot(data = data.frame(exprs.dat), genotypes = 
                              colData$genotypes)
    exprs.dat <- as.data.frame(exprs.dat)
    exprs.dat <- exprs.dat %>% rownames_to_column() %>% dplyr::rename(RefSeq.ID 
                                                                      = rowname)
    exprs.dat <- inner_join(x = exprs.dat, y = hg19.annot[,c("RefSeq.ID",
                        "gene.name","gene.length","transcript.length")], 
                        by = "RefSeq.ID")
    A.samples <- exprs.dat %>% dplyr::select(contains("SEQC_ILM_NVS_A")) %>% 
                                                                    colnames()
    B.samples <- exprs.dat %>% dplyr::select(contains("SEQC_ILM_NVS_B")) %>% 
                                                                    colnames()
    C.samples <- exprs.dat %>% dplyr::select(contains("SEQC_ILM_NVS_C")) %>% 
                                                                    colnames()
    dat.agg <- merge.dat.annot(exprs.dat)
    B.samples1 <- B.samples[c(1:32)]
    B.samples2 <- B.samples[c(33:64)]
    ratio.B1.B2 <- logofMeans.between.A.B(dat.agg, B.samples2, B.samples1)
    ratio.B1.B2 <- ratio.B1.B2[,c("logFC.crude","gene.length","gene.name")]
    Complot.B.seq <- gabels.plot(mat = data.frame(ratio.B1.B2), y.axis = 
                                   "Mean Log2 Fold Change")
    A.samples1 <- A.samples[c(1:16)]
    A.samples2 <- A.samples[c(49:64)]
    ratio.A1.A2 <- logofMeans.between.A.B(dat.agg, A.samples1, A.samples2)
    ratio.A1.A2 <- ratio.A1.A2[,c("logFC.crude","gene.length","gene.name")]
    Complot.A.seq <- gabels.plot(mat = data.frame(ratio.A1.A2), y.axis = 
                                   "Mean Log2 Fold Change")
    ratio.ABC.gl <- logofMeans.between.ABC(dat = dat.agg, A.samples, B.samples, 
                                           C.samples)
    ratio.ABC.gl <- data.frame(ratio.ABC.gl[, c("logFC.crude","gene.length", 
                                                "gene.name")])
    Complot.Beta.seq <- gabels.plot(mat = ratio.ABC.gl, y.axis = 
                                        "Mean Log2 (ß ratio)")
    
    
    featureData <- read.table(file = 
            "../dat-infofiles/annot.files/RefSeq.hg19.Gene.TransciptLength.txt", 
            sep = "\t", header = TRUE, quote = "", fill = TRUE,
            stringsAsFactors = FALSE, na.strings = FALSE)
    hg19.tx.annot <- aggregate(. ~ gene.name, data = 
                                   unique(featureData[,c("gene.name",
                                                         "transcript.length")]),
                               max)
    ratio.ABC.tx <- logofMeans.between.ABC(dat = dat.agg,A.samples,B.samples,
                                           C.samples)
    ratio.ABC.tx <- inner_join(x = ratio.ABC.tx, y = hg19.tx.annot,
                               by = "gene.name")
    Complot.Beta.seq.tx <- gabels.plot(mat=unique(ratio.ABC.tx[,c("logFC.crude",
                                                         "transcript.length")]), 
                                    length.type = "Transcript", 
                                    y.axis = "Mean Log2 (ß ratio)")
    seq.ABC.gl <- ratio.ABC.gl[,c("gene.name", "logFC.crude", "gene.length")]
    Complot.Beta.seq.tl <- gabels.plot(mat = ratio.ABC.gl, y.axis = 
                                        "Mean Log2 (ß ratio)")
    save(file = "../dat/SEQC/SEQC_RNA-Seq.RData", seq.ABC.gl)
    figures.array <- SEQC.microarray.analysis()
    comp.tech <- Complog2FC_tech()
    load("../results/SEQC_Brain-Rep_DEGs.RData")
    A.samples <- counts.table %>% dplyr::select(contains("SEQC_ILM_NVS_A")) %>% 
                                                                    colnames()
    B.samples <- counts.table %>% dplyr::select(contains("SEQC_ILM_NVS_B")) %>% 
                                                                    colnames()
    C.samples <- counts.table %>% dplyr::select(contains("SEQC_ILM_NVS_C")) %>% 
                                                                    colnames()
    res <- All.normalizationMethods(counts.table, A.samples, B.samples,
                                    C.samples, hg19.annot)
    B.samples1 <- B.samples[c(1:16)]
    B.samples2 <- B.samples[c(17:32)]
    B.samples3 <- B.samples[c(33:48)]
    B.samples4 <- B.samples[c(49:64)]
    samples.libID2 <- paste("Lib",samples.libID[c(1:48)],sep="")
    PCAplot.Brep <- PCAplot(data = exprs.dat[,c(B.samples1,B.samples2,B.samples3)],
                             genotypes = samples.libID2)
    ratio.AB.1 <- logofMeans.between.A.B(dat = dat.agg, B.samples = B.samples1, 
                                         A.samples = B.samples2)
    p1 <- gabels.plot(mat = ratio.AB.1[,c("logFC.crude", "gene.length")])
    ratio.AB.2 <- logofMeans.between.A.B(dat.agg, B.samples = B.samples3, 
                                         A.samples = B.samples2)
    p2 <- gabels.plot(mat = ratio.AB.2[,c("logFC.crude", "gene.length")])
    figures.seq <- list(figure4 = list(plotA = MDSplot.seq,
                                      plotB = Complot.Beta.seq,
                                      plotC = figures.array$figure4.array$plotC,
                                      plotD = figures.array$figure4.array$plotD,
                                      plotE = comp.tech$p5$plot,
                                      plotF = comp.tech$p6$plot,
                                      plotG = comp.tech$p7$plot),
                        figureS6.seq = list(plotA1 = Complot.B.seq,
                                            plotA2 = Complot.A.seq,
                                            plotA3 = Complot.Beta.seq.tx,
                                plotB1 = figures.array$figureS6.array$plotB1,
                                plotB2 = figures.array$figureS6.array$plotB2,
                                plotB3 = figures.array$figureS6.array$plotB3),
                        figureS7 = res, figureS8 = list(plotA = comp.tech$p1,
                                                        plotB = comp.tech$p2,
                                                        plotC = comp.tech$p3,
                                                        plotD = comp.tech$p4),
                        figureS9 = list(plotA = PCAplot.Brep, plotB = p1,
                                        plotC = p2, plotD1 = plots.supp9$p1,
                                        plotD2 = plots.supp9$p2,
                                        plotD3 = plots.supp9$p3,
                                        plotD4 = plots.supp9$p4))
    #save(file = "../results/SEQC_Figures.RData", figures.seq)
    return(figures.seq)
}