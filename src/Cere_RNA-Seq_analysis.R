################################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 19th April 2018
#
# Program is used for:
# 1. Diff Expression Analysis of Cerebellum RNA-Seq Dataset
################################################################################

cere_seq <- function(){

    ## files
    counts.table <- read.table(
                file = "../dat/Mecp2_Seq_NanoString/Cere.counts.table_mm10.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    rownames(counts.table) <- counts.table$ENSEMBL.ID
    counts.table$gene.length <- counts.table$end - counts.table$start + 1
    
    ## coldData and dataset with all the gene Ids
    genotypes <- factor(c(rep("KO", 3), rep("WT", 3)), levels = c("KO", "WT"))
    colData <- data.frame(row.names = colnames(counts.table[, c(2:7)]), 
                          genotypes = genotypes)
    head(colData)
    dim(colData)
    
    ## DESeq2
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts.table[, c(2:7)]),
                                  colData = colData, design = ~ genotypes)
    dds <- estimateSizeFactors(dds)
    dat <- counts(dds, normalized = TRUE)
    
    ## Box and PCA plot for normalized dataset
    #boxPlot(data = log2(dat + 1), title = "", samples = colData$genotypes)
    p1 <- PCAplot(data = log2(dat + 1), genotypes)
    
    ## results from DESeq2
    dds <- DESeq(dds, betaPrior = TRUE)
    res.KO <- results(dds, contrast = c("genotypes", "KO", "WT"))
    res.KO$norm.counts <- counts(dds, normalized = TRUE)
    res.KO.annot <- inner_join(x = data.frame(res.KO) %>% 
                                    rownames_to_column(var = "ENSEMBL.ID"),
                                y = counts.table[, c(1, 8:14)],
                                by = "ENSEMBL.ID")
    
    ## Scatter Plot
    res.KO.annot <- logofMeans.between.A.B(res.KO.annot, A.samples = c(11:13),
                                           B.samples = c(8:10))
    p2 <- Plot.Scatter(dat = res.KO.annot[, c("gene.name", "log2FoldChange", 
                                              "padj", "gene.length")],
                       log2FC = log2(1.15), comp.between =  "(KO/WT) RNA-Seq")
    
    ## Overlay plots
    res.KO.annot$comp.mat1 <- apply(X = data.frame(res.KO.annot[,c(11:13)]), 1,
                                    function(r) {log2((r[3] + 1)/(r[1] + 1))})
    p3 <- overlay.gabels.plot(mat = res.KO.annot[,c("comp.mat1", "logFC.crude",
                                                    "gene.length")],
                              comp.between1 = "(WT/WT)",
                              comp.between2 = "(KO/WT)")
    p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.25, 0.25)),
                    p3$plot2 + coord_cartesian(ylim = c(0, 20)), ncol = 1,
                    align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5),
                                                            "cm"))
    
    # ## writing whole genome DEGs
    # write.table(x = data.frame(res.KO) %>% rownames_to_column(var = 
    #                                                            "ENSEMBL.ID"),
    #       file = "../dat/Mecp2_Seq_NanoString/Cere_Seq_WholeGenome_DEGs.txt",
    #       quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    # write.table(x = res.KO.annot, file = 
    #               "../dat/Mecp2_Seq_NanoString/Cere_Seq_WholeGenome_Comp.txt",
    #               quote = FALSE, sep="\t", row.names=FALSE,col.names = TRUE)
    
    ## genes present in Nanostring only
    Mecp2.ns <- read.table("../dat/Mecp2_Seq_NanoString/NS_Mecp2.txt",
                           sep = "\t", stringsAsFactors = FALSE, quote = "",
                           header = TRUE, na.strings = FALSE)
    gene.ns <- which(Mecp2.ns$Code.Class == "Endogenous")
    gene.ns <- Mecp2.ns[gene.ns, "Name"]
    
    ## RNA-Seq with NS gene only
    counts.table[which(counts.table$gene.name == "9430076C15Rik"), 
                 "gene.name"] <- "Creb5"
    counts.table <- counts.table[which(counts.table$gene.name %in% gene.ns), ]
    rownames(counts.table) <- counts.table$gene.name
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts.table[, c(2:7)]),
                                  colData = colData, design = ~ genotypes)
    dds <- estimateSizeFactors(dds)
    dat <- counts(dds, normalized = TRUE)
    p4 <- PCAplot(data = log2(dat + 1), genotypes)
    dds <- DESeq(dds, betaPrior = TRUE)
    res.KO <- results(dds, contrast = c("genotypes", "KO", "WT"))
    res.KO$norm.counts <- counts(dds, normalized = TRUE)
    res.KO <- data.frame(res.KO)
    res.KO$logFC.crude <- apply(res.KO[,c(7:12)],1,function(r)
                                {log2((mean(r[1:3])+1)/(mean(r[4:6])+1))})
    res.KO.annot <- inner_join(x = res.KO[,c(2,13,6)] %>% 
                                   rownames_to_column(var = "gene.name"),
                               y = counts.table[,c(8:14)],
                               by = "gene.name")
    res.KO.annot <- res.KO.annot[order(res.KO.annot$gene.name, decreasing = F),]
    plots.cere.seq <- list(p1 = p1, p2 = p3, p3 = p2, p4 = p4, seq.annot = 
                               res.KO.annot)
    return(plots.cere.seq)
}