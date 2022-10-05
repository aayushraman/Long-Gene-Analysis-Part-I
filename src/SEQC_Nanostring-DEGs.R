#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 23rd April 2017
#
# Program is used for:
# 1. SEQC dataset
# 2. Comparison with SEQC Nanostring 
#############################################################################

rm(list = ls())

## functions and set the working directory
setwd("~/Desktop/AayushRaman/RWorkspace/MeCP2-HudaZoghbi/LongGene-ShortGeneAnalysis/code_manuscript/src_version0.1/")
source("SEQC_libraries.R")

## SEQC RNA-Seq Dataset
counts.table <- read.table(file = "../dat/SEQC/GSE47774_SEQC_ILM_NVS.txt", header = TRUE, sep = "\t", 
                           quote = "", row.names = 1, stringsAsFactors = FALSE)
hg19.annot <- read.table(file = "../dat-infofiles/annot.files/RefSeq.hg19.Gene.TransciptLength.txt", header = TRUE, sep = "\t", 
                         quote = "", stringsAsFactors = FALSE)
erccindex <- grep("ERCC",rownames(counts.table))
counts.table <- counts.table[-erccindex,]

## Nanostring
seqc.NanoString <- read.table("../dat/SEQC/nCounterDataset.txt", sep = "\t", 
                              stringsAsFactors = FALSE, quote = "", header = TRUE, na.strings = FALSE)
gene.ns <- seqc.NanoString$Name[which(seqc.NanoString$Code.Class == "Endogenous")]    

## sample types
samples.type <- gsub(pattern = "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_(.*)_(.*)", 
                     replacement = "\\3", x = colnames(counts.table))
samples.barcode <- gsub(pattern = "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_([ATGC]+)_(.*)", 
                        replacement = "\\6", x = colnames(counts.table))
samples.libID <- gsub(pattern = "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_([ATGC]+)_(.*)", 
                      replacement = "\\4", x = colnames(counts.table))
samples.flowcell <- gsub(pattern = "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_([ATGC]+)_(.*)", 
                         replacement = "\\7", x = colnames(counts.table))

## merging the annotation file
counts.table <- counts.table %>% rownames_to_column(var = "RefSeq.ID")
counts.table <- inner_join(x = counts.table, y = hg19.annot, by = "RefSeq.ID")
counts.table <- counts.table[which(counts.table$gene.name %in% gene.ns),]
rownames(counts.table) <- counts.table$RefSeq.ID
counts.table <- counts.table[,-1]

## Arguments for DESeq2 and its run
## Including the genes present in Nanostring only
ind <- which(samples.type == "A" | samples.type == "B" | samples.type == "C" | samples.type == "D")
counts.table <- counts.table[,ind]
colData <- as.data.frame(samples.type[ind])
rownames(colData) <- colnames(counts.table)
colnames(colData) <- "genotypes"

##############################
#
## DESeq2 on SEQC dataset
#
##############################

## Tg/WT
idx <- which(colData$genotypes == "A" | colData$genotypes == "D")
dat <- counts.table[,idx]
colData1 <- data.frame(row.names = colnames(dat), 
                        genotypes = colData[idx,])
dds.seqc <- DESeqDataSetFromMatrix(countData = as.matrix(dat), colData1, design = ~ genotypes)
dds.seqc <- DESeq(dds.seqc)
res.D.A <- results(dds.seqc, contrast = c("genotypes", "A", "D"))

## KO/WT
idx <- which(colData$genotypes == "B" | colData$genotypes == "D")
dat <- counts.table[,idx]
colData1 <- data.frame(row.names = colnames(dat), 
                       genotypes = colData[idx,])
dds.seqc <- DESeqDataSetFromMatrix(countData = as.matrix(dat), colData1, design = ~ genotypes)
dds.seqc <- DESeq(dds.seqc)
res.D.B <- results(dds.seqc, contrast = c("genotypes", "B", "D"))

## inner_join with hg19.annot
res.D.A$RefSeq.ID <- rownames(res.D.A)
res.D.A$logFC.crude <- res.D.A$log2FoldChange
res.D.B$RefSeq.ID <- rownames(res.D.B)
res.D.B$logFC.crude <- res.D.B$log2FoldChange
res.D.A <- inner_join(x = data.frame(res.D.A), y = hg19.annot, by = "RefSeq.ID") 
res.D.B <- inner_join(x = data.frame(res.D.B), y = hg19.annot, by = "RefSeq.ID")

## remove the duplicated genes
dup.degs.1 <- unique(res.D.A[duplicated(res.D.A$gene.name),"gene.name"])
dup.degs.2 <- unique(res.D.B[duplicated(res.D.B$gene.name),"gene.name"])
nmID.rm <- c()
remove.dups <- function(mat, dups){
    for(i in 1:length(dups)){
        gene <- dups[i]
        idx <- which(mat$gene.name == gene)
        rm.idx <- idx[which(!idx %in% idx[which.max(abs(mat$log2FoldChange)[idx])])]
        nmID.rm <- c(nmID.rm, rm.idx)
    }
    mat <- mat[-nmID.rm, ]
}
res.D.A <- remove.dups(res.D.A, dup.degs.1)
res.D.B <- remove.dups(res.D.B, dup.degs.1)

## Scatter plot lm
p1 <- scatter.lm(unique(res.D.A[,c("gene.name", "logFC.crude", "gene.length")])) 
p1 + coord_cartesian(ylim = c(-5,2)) + geom_hline(yintercept = 0)

p2 <- scatter.lm(unique(res.D.B[,c("gene.name", "logFC.crude", "gene.length")])) 
p2 + coord_cartesian(ylim = c(-5,2)) + geom_hline(yintercept = 0)

## DEGs Scatter plot 
Plot.Scatter(dat = res.D.A[,c("gene.name", "log2FoldChange", "padj", "gene.length")], 
             log2FC = log2(1.15), comp.between =  "(Tg/WT)")
Plot.Scatter(dat = res.D.B[,c("gene.name", "log2FoldChange", "padj", "gene.length")], 
             log2FC = log2(1.15), comp.between =  "(KO/WT)")

## Nanostring
seqc.NanoString.norm <- NanoStringNorm(x = seqc.NanoString, anno = NA,CodeCount = 'sum', Background = 'mean',
                                           SampleContent = 'housekeeping.sum', round.values = FALSE, take.log = FALSE,
                                           return.matrix.of.endogenous.probes = TRUE)
genotypes.ns <- factor(rep(c(rep("A",3),rep("B",3),rep("C",3),rep("D",3)),2), levels = c("A", "B", "C","D"))
idx <- order(genotypes.ns)
genotypes.ns <- genotypes.ns[idx]
seqc.NanoString.norm <- seqc.NanoString.norm[,idx]

## PCA Plot
cat("\n\n Printing Supp Figure 3(A)\n\n")
boxPlot(data = log2(seqc.NanoString.norm+1), title = "", samples = genotypes.ns)    
MDSplot(data = log2(seqc.NanoString.norm+1), genotypes = genotypes.ns)

## hg19
hg19.ns <- hg19.annot[,c("gene.name", "gene.length")]
hg19.ns <- hg19.ns[which(hg19.ns$gene.name %in% rownames(seqc.NanoString.norm)),]
hg19.ns <- aggregate(gene.length ~ gene.name, data = hg19.ns, max)

## DEGs analysis -- KO
seqc.NanoString.norm <- data.frame(seqc.NanoString.norm)
dds.KO <- DESeqDataSetFromMatrix(countData = as.matrix(round(seqc.NanoString.norm[,c(7:12,19:24)])),
                                 colData = data.frame(row.names = colnames(seqc.NanoString.norm[,c(7:12,19:24)]),
                                                      genotypes = genotypes.ns[which(genotypes.ns %in% c("B","D"))]),
                                 design = ~ genotypes)
sizeFactors(dds.KO) <- rep(1,12)
dds.KO <- DESeq(dds.KO)
res.KO <- results(dds.KO, contrast = c("genotypes", "B", "D"))
res.KO$gene <- rownames(res.KO)
res.KO <- inner_join(x = data.frame(res.KO), y = hg19.ns, by = c("gene" = "gene.name"))
sum(res.KO$padj < 0.05, na.rm = TRUE) ## 502

## DEGs analysis -- Tg3
dds.Tg <- DESeqDataSetFromMatrix(countData = as.matrix(round(seqc.NanoString.norm[,c(1:6,19:24)])),
                                 colData = data.frame(row.names = colnames(seqc.NanoString.norm[,c(1:6,19:24)]),
                                                      genotypes = genotypes.ns[which(genotypes.ns %in% c("A","D"))]),
                                 design = ~ genotypes)
sizeFactors(dds.Tg) <- rep(1,12)
dds.Tg <- DESeq(dds.Tg)
res.Tg <- results(dds.Tg, contrast = c("genotypes", "A", "D"))
res.Tg$gene <- rownames(res.Tg)
res.Tg <- inner_join(x = data.frame(res.Tg), y = hg19.ns, by = c("gene" = "gene.name"))
sum(res.Tg$padj < 0.05, na.rm = TRUE) ## 610

## DEGs Scatter plot
logFC.comp <- function(idx, fc, RNASeq.logFC, NS.logFC, RNASeq.fdr, NS.fdr){
    res.KO.HiC <- res.Tg.HiC 
    dat <- res.Tg.HiC[idx,] #res.KO.HiC[idx,]
    mat1 <- dat[which(dat$longGene == "LongGene"),]
    mat2 <- dat[which(dat$longGene == "ShortGene"),]
    q1 <- qplot(data = mat1, x = log2FoldChange.Seq, y = log2FoldChange.NS, col = Col.degs) + 
        geom_point(shape=1) + geom_jitter() + geom_abline(colour = "grey", size = 0.5) + 
        xlim(c(-4, 2)) + ylim(c(-4, 2)) + ggtitle("Long Genes") + 
        geom_hline(aes(yintercept=0), linetype="dotted") + 
        geom_vline(aes(xintercept=0), linetype="dotted") +        
        geom_hline(aes(yintercept=log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=log2(1.2)), linetype="dotted") +
        geom_hline(aes(yintercept=-log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=-log2(1.2)), linetype="dotted") +
        xlab(paste("LogFC RNA-Seq Values")) + ylab(paste("LogFC NanoString Values")) + 
        theme(legend.text = element_text(size = 12, face = "bold"), 
              legend.title = element_text(size = 12, colour = "black", face = "bold"),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(size = 18, face = "bold", color = "black"),
              axis.text.y = element_text(size = 18, face = "bold", color = "black"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    print(q1)
    
    q2 <- qplot(data = mat2, x = log2FoldChange.Seq, y = log2FoldChange.NS, col = Col.degs) + 
        geom_point(shape=1) + geom_jitter() + geom_abline(colour = "grey", size = 0.5) + 
        xlim(c(-4, 2)) + ylim(c(-4, 2)) + ggtitle("Short Genes") + 
        geom_hline(aes(yintercept=0), linetype="dotted") + 
        geom_vline(aes(xintercept=0), linetype="dotted") +        
        geom_hline(aes(yintercept=log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=log2(1.2)), linetype="dotted") +
        geom_hline(aes(yintercept=-log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=-log2(1.2)), linetype="dotted") +
        xlab(paste("LogFC RNA-Seq Values")) + ylab(paste("LogFC NanoString Values")) + 
        theme(legend.text = element_text(size = 12, face = "bold"), 
              legend.title = element_text(size = 12, colour = "black", face = "bold"),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(size = 18, face = "bold", color = "black"),
              axis.text.y = element_text(size = 18, face = "bold", color = "black"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    print(q2)
    
    ### HiC Dataset 
    cat("RNA-Seq Dataset\n")
    cat("Total Number of DEGs = ", sum(RNASeq.fdr < 0.05 & abs(RNASeq.logFC) > log2(fc), na.rm = TRUE),"\n")
    cat("Total Number of Up.DEGs = ", sum(RNASeq.fdr < 0.05 & RNASeq.logFC > log2(fc), na.rm = TRUE),"\n")
    cat("Total Number of Down.DEGs = ", sum(RNASeq.fdr < 0.05 & RNASeq.logFC < -log2(fc), na.rm = TRUE),"\n")
    
    ### HiC Long Gene Dataset 
    cat("Total Number of Long DEGs = ", sum(RNASeq.fdr < 0.05 & abs(RNASeq.logFC) > log2(fc) &
                                                res.KO.HiC$longGene == "LongGene", na.rm = TRUE),"\n")
    cat("Total Number of Long Up DEGs = ", sum(RNASeq.fdr < 0.05 & RNASeq.logFC > log2(fc) &
                                                   res.KO.HiC$longGene == "LongGene", na.rm = TRUE),"\n")
    cat("Total Number of Long Down DEGs = ", sum(RNASeq.fdr < 0.05 & RNASeq.logFC < -log2(fc) &
                                                     res.KO.HiC$longGene == "LongGene", na.rm = TRUE),"\n")
    
    ### HiC Short Gene Dataset 
    cat("Total Number of Short DEGs = ", sum(RNASeq.fdr < 0.05 & abs(RNASeq.logFC) > log2(fc) &
                                                 res.KO.HiC$longGene == "ShortGene", na.rm = TRUE),"\n")
    cat("Total Number of Short Up DEGs = ", sum(RNASeq.fdr < 0.05 & RNASeq.logFC > log2(fc) &
                                                    res.KO.HiC$longGene == "ShortGene", na.rm = TRUE),"\n")
    cat("Total Number of Short Down DEGs = ", sum(RNASeq.fdr < 0.05 & RNASeq.logFC < -log2(fc) &
                                                      res.KO.HiC$longGene == "ShortGene", na.rm = TRUE),"\n")
    
    ### NanoString Dataset
    cat("NanoString Dataset\n\n")
    cat("Total Number of DEGs = ", sum(NS.fdr < 0.05 & abs(NS.logFC) > log2(fc), na.rm = TRUE),"\n")
    cat("Total Number of Up DEGs = ", sum(NS.fdr < 0.05 & NS.logFC > log2(fc), na.rm = TRUE),"\n")
    cat("Total Number of Down DEGs = ", sum(NS.fdr < 0.05 & NS.logFC < -log2(fc), na.rm = TRUE),"\n")
    
    cat("Total Number of Long DEGs = ", sum(NS.fdr < 0.05 & abs(NS.logFC) > log2(fc) &
                                                res.KO.HiC$longGene == "LongGene", na.rm = TRUE),"\n")
    cat("Total Number of Long Up DEGs = ", sum(NS.fdr < 0.05 & NS.logFC > log2(fc) &
                                                   res.KO.HiC$longGene == "LongGene", na.rm = TRUE),"\n")
    cat("Total Number of Long Down DEGs = ", sum(NS.fdr < 0.05 & NS.logFC < -log2(fc) &
                                                     res.KO.HiC$longGene == "LongGene", na.rm = TRUE),"\n")
    
    cat("Total Number of Short DEGs = ", sum(NS.fdr < 0.05 & abs(NS.logFC) > log2(fc) &
                                                 res.KO.HiC$longGene == "ShortGene", na.rm = TRUE),"\n")
    cat("Total Number of Short Up DEGs = ", sum(NS.fdr < 0.05 & NS.logFC > log2(fc) &
                                                    res.KO.HiC$longGene == "ShortGene", na.rm = TRUE),"\n")
    cat("Total Number of Short Down DEGs = ", sum(NS.fdr < 0.05 & NS.logFC < -log2(fc) &
                                                      res.KO.HiC$longGene == "ShortGene", na.rm = TRUE),"\n")
    
    ## Overlap
    cat("Total Number of Overlap Long Up DEGs = ", sum(res.KO.HiC$Col.degs == "Both" & res.KO.HiC$longGene == "LongGene" &
                                                           NS.logFC > log2(fc) & RNASeq.logFC > log2(fc), na.rm = TRUE),"\n")
    cat("Total Number of Overlap Long Down DEGs = ", sum(res.KO.HiC$Col.degs == "Both" & res.KO.HiC$longGene == "LongGene" &
                                                             NS.logFC < -log2(fc) & RNASeq.logFC < -log2(fc), na.rm = TRUE),"\n")
    
    cat("Total Number of Overlap Short Up DEGs = ", sum(res.KO.HiC$Col.degs == "Both" & res.KO.HiC$longGene == "ShortGene" &
                                                            NS.logFC > log2(fc) & RNASeq.logFC > log2(fc), na.rm = TRUE),"\n")
    cat("Total Number of Overlap Short down DEGs = ", sum(res.KO.HiC$Col.degs == "Both" & res.KO.HiC$longGene == "ShortGene" &
                                                              NS.logFC < -log2(fc) & RNASeq.logFC < -log2(fc), na.rm = TRUE),"\n")
}

## Scatter plot
res.KO$logFC.crude <- res.KO$log2FoldChange
res.D.B$logFC.crude <- res.D.B$log2FoldChange

res.Tg$logFC.crude <- res.Tg$log2FoldChange
res.D.A$logFC.crude <- res.D.A$log2FoldChange

## scatter plot
p1 <- scatter.lm(unique(res.KO[,c("gene", "logFC.crude", "gene.length")])) 
p1$plot
p1$plot + coord_cartesian(ylim = c(-5,2)) + geom_hline(yintercept = 0)

p2 <- scatter.lm(unique(res.Tg[,c("gene", "logFC.crude", "gene.length")])) 
p2$plot
p2$plot + coord_cartesian(ylim = c(-5,2)) + geom_hline(yintercept = 0)

## B/D plots
res.KO.HiC <- inner_join(x = res.D.B, y = res.KO, by = c("gene.name" = "gene"))
res.KO.HiC$longGene <- ifelse(res.KO.HiC$gene.length.x > 100e3, "LongGene", "ShortGene")
res.KO.HiC <- res.KO.HiC[,c("gene.name", "log2FoldChange.x", "log2FoldChange.y", "padj.x", "padj.y", "longGene", "gene.length.x")]
colnames(res.KO.HiC) <- c("gene.name", "log2FoldChange.Seq", "log2FoldChange.NS", "padj.Seq", 
                            "padj.NS", "longGene", "gene.length")
res.KO.HiC$Col.degs <- ifelse(res.KO.HiC$padj.Seq < 0.05 & abs(res.KO.HiC$log2FoldChange.Seq) > 0 &
                                  res.KO.HiC$padj.NS < 0.05 & abs(res.KO.HiC$log2FoldChange.NS) > 0,"Both",
                              ifelse(res.KO.HiC$padj.Seq < 0.05 & abs(res.KO.HiC$log2FoldChange.Seq) > 0, "RNA-Seq",
                                     ifelse(res.KO.HiC$padj.NS < 0.05 & abs(res.KO.HiC$log2FoldChange.Seq) > 0, "NanoString","Not Sig")))
sum(res.KO.HiC$longGene == "LongGene")
sum(res.KO.HiC$longGene == "ShortGene")
sum(res.KO.HiC$Col.degs == "RNA-Seq", na.rm = TRUE)
sum(res.KO.HiC$Col.degs == "NanoString", na.rm = TRUE)
sum(res.KO.HiC$Col.degs == "Both", na.rm = TRUE)

## AR_DESeq2 logFC -- KO
idx <- which(!is.na(res.KO.HiC$Col.degs) & res.KO.HiC$Col.degs != "Not Sig")
logFC.comp(idx = idx, fc = 1, RNASeq.logFC = res.KO.HiC$log2FoldChange.Seq, 
           NS.logFC = res.KO.HiC$log2FoldChange.NS, RNASeq.fdr = res.KO.HiC$padj.Seq, NS.fdr = res.KO.HiC$padj.NS)
logFC.comp(idx = idx, fc = 1.1, RNASeq.logFC = res.KO.HiC$log2FoldChange.Seq, 
           NS.logFC = res.KO.HiC$log2FoldChange.NS, RNASeq.fdr = res.KO.HiC$padj.Seq, NS.fdr = res.KO.HiC$padj.NS)
logFC.comp(idx = idx, fc = 1.15, RNASeq.logFC = res.KO.HiC$log2FoldChange.Seq, 
           NS.logFC = res.KO.HiC$log2FoldChange.NS, RNASeq.fdr = res.KO.HiC$padj.Seq, NS.fdr = res.KO.HiC$padj.NS)
logFC.comp(idx = idx, fc = 1.2, RNASeq.logFC = res.KO.HiC$log2FoldChange.Seq, 
           NS.logFC = res.KO.HiC$log2FoldChange.NS, RNASeq.fdr = res.KO.HiC$padj.Seq, NS.fdr = res.KO.HiC$padj.NS)

Plot.Scatter(dat = res.KO.HiC[,c("gene.name", "log2FoldChange.NS", "padj.NS", "gene.length")], 
                log2FC = log2(1.15), comp.between =  "(KO/WT)")
rm(res.KO.HiC)

## A/D plots
res.Tg.HiC <- inner_join(x = res.D.A, y = res.Tg, by = c("gene.name" = "gene"))
res.Tg.HiC$longGene <- ifelse(res.Tg.HiC$gene.length.x > 100e3, "LongGene", "ShortGene")
res.Tg.HiC <- res.Tg.HiC[,c("gene.name", "log2FoldChange.x", "log2FoldChange.y", "padj.x", "padj.y", "longGene", "gene.length.x")]
colnames(res.Tg.HiC) <- c("gene.name", "log2FoldChange.Seq", "log2FoldChange.NS", "padj.Seq", 
                          "padj.NS", "longGene", "gene.length")
res.Tg.HiC$Col.degs <- ifelse(res.Tg.HiC$padj.Seq < 0.05 & abs(res.Tg.HiC$log2FoldChange.Seq) > 0 &
                                  res.Tg.HiC$padj.NS < 0.05 & abs(res.Tg.HiC$log2FoldChange.NS) > 0,"Both",
                              ifelse(res.Tg.HiC$padj.Seq < 0.05 & abs(res.Tg.HiC$log2FoldChange.Seq) > 0, "RNA-Seq",
                                     ifelse(res.Tg.HiC$padj.NS < 0.05 & abs(res.Tg.HiC$log2FoldChange.Seq) > 0, "NanoString","Not Sig")))
sum(res.Tg.HiC$longGene == "LongGene")
sum(res.Tg.HiC$longGene == "ShortGene")
sum(res.Tg.HiC$Col.degs == "RNA-Seq", na.rm = TRUE)
sum(res.Tg.HiC$Col.degs == "NanoString", na.rm = TRUE)
sum(res.Tg.HiC$Col.degs == "Both", na.rm = TRUE)

## AR_DESeq2 logFC -- Tg
idx <- which(!is.na(res.Tg.HiC$Col.degs) & res.Tg.HiC$Col.degs != "Not Sig")
logFC.comp(idx = idx, fc = 1, RNASeq.logFC = res.Tg.HiC$log2FoldChange.Seq, 
           NS.logFC = res.Tg.HiC$log2FoldChange.NS, RNASeq.fdr = res.Tg.HiC$padj.Seq, NS.fdr = res.Tg.HiC$padj.NS)
logFC.comp(idx = idx, fc = 1.1, RNASeq.logFC = res.Tg.HiC$log2FoldChange.Seq, 
           NS.logFC = res.Tg.HiC$log2FoldChange.NS, RNASeq.fdr = res.Tg.HiC$padj.Seq, NS.fdr = res.Tg.HiC$padj.NS)
logFC.comp(idx = idx, fc = 1.15, RNASeq.logFC = res.Tg.HiC$log2FoldChange.Seq, 
           NS.logFC = res.Tg.HiC$log2FoldChange.NS, RNASeq.fdr = res.Tg.HiC$padj.Seq, NS.fdr = res.Tg.HiC$padj.NS)
logFC.comp(idx = idx, fc = 1.2, RNASeq.logFC = res.Tg.HiC$log2FoldChange.Seq, 
           NS.logFC = res.Tg.HiC$log2FoldChange.NS, RNASeq.fdr = res.Tg.HiC$padj.Seq, NS.fdr = res.Tg.HiC$padj.NS)

Plot.Scatter(dat = res.Tg.HiC[,c("gene.name", "log2FoldChange.NS", "padj.NS", "gene.length")], 
             log2FC = log2(1.15), comp.between =  "(Tg/WT)")
rm(res.KO.HiC)

