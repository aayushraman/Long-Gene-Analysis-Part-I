#############################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 18th April 2018
#
# Program is used for:
# 1. DEGs on Brain Samples based on technical replicates
#############################################################################

rm(list = ls())

## functions and set the working directory
source("libraries.R")

## NVS dataset and hg19 annot.
ns.genes <- read.table(file = "../dat/SEQC/PanCanHumanGenes.txt", header = FALSE, 
                       sep = "\t", quote = "", stringsAsFactors = FALSE, 
                       col.names = c("gene"))
counts.table <- read.table(file = "../dat/SEQC/GSE47774_SEQC_ILM_NVS.txt", 
                           header = TRUE, sep = "\t",quote = "", row.names = 1, 
                           stringsAsFactors = FALSE)
hg19.annot <- read.table(file = "../dat-infofiles/annot.files/RefSeq.hg19.Gene.TransciptLength.txt",
                         header = TRUE, sep = "\t", quote = "", 
                         stringsAsFactors = FALSE)
erccindex <- grep("ERCC",rownames(counts.table))
counts.table <- counts.table[-erccindex,]
  
## sample types
samples.type <- gsub(pattern = "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_(.*)_(.*)", 
                       replacement = "\\3", x = colnames(counts.table))
samples.libID <- gsub(pattern = "SEQC_(ILM|LIF)_(.*)_([A-Z])_([0-9]+)_L([0-9]+)_([ATGC]+)_(.*)", 
                        replacement = "\\4", x = colnames(counts.table))

## B Samples
B.samples <- counts.table %>% dplyr::select(contains("SEQC_ILM_NVS_B")) %>% colnames()
KO <- B.samples[c(1:16)]
WT <- B.samples[c(17:32)]
Tg <- B.samples[c(33:48)]
counts.table <- counts.table[,c(KO,WT,Tg)]
colData <- data.frame(row.names = colnames(counts.table), 
                      genotypes = c(rep("KO",16), rep("WT",16),rep("Tg3",16)))

## dds human genes dataset
exprs.dat <- DESeqCalculation(seqcData = counts.table, genotype = colData)
MDSplot(data = data.frame(exprs.dat), genotypes = colData$genotypes)
PCAplot(data = data.frame(exprs.dat), genotypes = colData$genotypes)

## long gene bias plots
exprs.dat <- as.data.frame(exprs.dat)
exprs.dat <- exprs.dat %>% rownames_to_column() %>% dplyr::rename(RefSeq.ID = rowname)
exprs.dat <- inner_join(x = exprs.dat, y = hg19.annot[,c("RefSeq.ID","gene.name",
                        "gene.length","transcript.length")], by = "RefSeq.ID")
dat.agg <- merge.dat.annot(exprs.dat)
ratio.AB.1 <- logofMeans.between.A.B(dat = dat.agg, B.samples = KO, A.samples = WT)
ratio.AB.2 <- logofMeans.between.A.B(dat.agg, B.samples = Tg, A.samples = WT)
p1 <- gabels.plot(mat = ratio.AB.1[,c("logFC.crude", "gene.length")])
p2 <- gabels.plot(mat = ratio.AB.2[,c("logFC.crude", "gene.length")])

## DEGs Analysis -- KO
counts.KO <- counts.table[,c(KO,WT)] %>% rownames_to_column() %>% 
                dplyr::rename(RefSeq.ID = rowname)
counts.KO <- inner_join(x = counts.KO, y = hg19.annot[,c("RefSeq.ID",
            "gene.name","gene.length","transcript.length")], by = "RefSeq.ID")
counts.KO <- merge.dat.annot(counts.KO)
rownames(counts.KO) <- counts.KO$gene.name
coldata.KO <- data.frame(row.names = colnames(counts.KO[2:33]), 
                            genotypes = c(rep("KO",16), rep("WT",16)))
dds.KO <- DESeqDataSetFromMatrix(countData = counts.KO[which(counts.KO$gene.name 
                                                             %in% ns.genes$gene),c(2:33)], 
                                 colData = coldata.KO, design = ~ genotypes)
dds.KO <- DESeq(dds.KO)
res.KO <- results(dds.KO, contrast = c("genotypes", "KO", "WT"))
res.KO$gene.name <- rownames(res.KO)
res.KO <- inner_join(x = data.frame(res.KO), y = counts.KO[,c("gene.name",
                                                              "gene.length")], 
                     by = "gene.name")
res.KO$longGene <- ifelse(res.KO$gene.length > 100e3, "LongGene", "ShortGene")

## DEGs Scatter Plot 
Plot.Scatter <- function(dat, log2FC, comp.between, pval = 0.05){
    colnames(dat) <- c("gene.name", "logFC", "adj.P.Val", "gene.length")
    gene.type <- ifelse((dat$adj.P.Val < pval & abs(dat$logFC) > log2FC), 
                        ifelse(dat$gene.length > 100e3, "Long Genes", "Short Genes"),"Not Stat. Signif.")
    ind <- which(gene.type != "Not Stat. Signif.")
    dat <- dat[ind,]
    gene.type <- gene.type[ind]
    
    ## Contingency Table
    up.LongGenes <- sum(dat$logFC > log2FC & dat$gene.length > 100e3)
    down.LongGenes <- sum(dat$logFC < -log2FC & dat$gene.length > 100e3)
    up.ShortGenes <- sum(dat$logFC > log2FC & dat$gene.length <= 100e3)
    down.ShortGenes <- sum(dat$logFC < -log2FC & dat$gene.length <= 100e3)
    cont.tab <- matrix(data = c(up.LongGenes, down.LongGenes, up.ShortGenes, down.ShortGenes), nrow = 2, ncol = 2)
    rownames(cont.tab) <- c("Up", "Down")
    colnames(cont.tab) <- c("Long.Gene", "Short.Gene")
    print(cont.tab)
    
    ## qplot
    print(qplot(y = dat$logFC, x = dat$gene.length/1000, colour = gene.type,
                xlab = "Gene Length in KB", ylab = paste("Log2 Fold Change", comp.between)) + 
              scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) + 
              coord_cartesian(ylim = c(-1.5,1.5)) +
              annotate("text", x = 250, y=0.25, label= cont.tab[1,1], size=7, fontface="bold") + 
              annotate("text", x = 250, y=-0.25, label= cont.tab[2,1], size=7, fontface="bold") + 
              annotate("text", x = 5, y=0.25, label = cont.tab[1,2], size=7, fontface="bold") +
              annotate("text", x = 5, y=-0.25, label = cont.tab[2,2], size=7, fontface="bold") + theme_grey()+
              theme(plot.title = element_text(size = 14, face = "bold"),
                    axis.title = element_text(size = 18, face = "bold"),
                    legend.position="none", axis.text.x = element_text(size = 18, face = "bold", color = "black"),
                    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
                    legend.text = element_text(size = 18, face = "bold"),
                    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
}
p1 <- Plot.Scatter(dat = res.KO[,c("gene.name", "log2FoldChange", "padj", "gene.length")], 
             log2FC = log2(1), comp.between =  "(Lib1/Lib2)")
p1 <- p1 + coord_cartesian(ylim = c(-0.3,0.3))

p2 <- Plot.Scatter(dat = res.KO[,c("gene.name", "log2FoldChange", "padj", "gene.length")], 
                   log2FC = log2(1.05), comp.between =  "(Lib1/Lib2)")
p2 <- p2 + coord_cartesian(ylim = c(-0.3,0.3))

## DEGs Analysis -- Tg
counts.Tg <- counts.table[,c(Tg,WT)] %>% rownames_to_column() %>% 
                dplyr::rename(RefSeq.ID = rowname)
counts.Tg <- inner_join(x = counts.Tg, y = hg19.annot[,c("RefSeq.ID","gene.name",
                            "gene.length","transcript.length")], by = "RefSeq.ID")
counts.Tg <- merge.dat.annot(counts.Tg)
rownames(counts.Tg) <- counts.Tg$gene.name
coldata.Tg <- data.frame(row.names = colnames(counts.Tg[2:33]), genotypes = c(rep("Tg",16), rep("WT",16)))
dds.Tg <- DESeqDataSetFromMatrix(countData = counts.Tg[which(counts.Tg$gene.name %in% ns.genes$gene),
                                                       c(2:33)], colData = coldata.Tg, design = ~ genotypes)
dds.Tg <- DESeq(dds.Tg)
res.Tg <- results(dds.Tg, contrast = c("genotypes", "Tg", "WT"))
res.Tg$gene.name <- rownames(dds.Tg)
sum(res.Tg$padj < 0.05, na.rm = TRUE)
sum(res.Tg$log2FoldChange > log2(1.15) & dds.Tg$padj < 0.05, na.rm = TRUE)
sum(res.Tg$log2FoldChange < -log2(1.15) & dds.Tg$padj < 0.05, na.rm = TRUE)
res.Tg <- inner_join(x = data.frame(res.Tg), y = counts.KO[,c("gene.name","gene.length")], by = "gene.name")
res.Tg$longGene <- ifelse(res.KO$gene.length > 100e3, "LongGene", "ShortGene")

## DEGs Scatter Plot 
p3 <- Plot.Scatter(dat = res.Tg[,c("gene.name", "log2FoldChange", "padj", "gene.length")], 
                   log2FC = log2(1), comp.between =  "(Lib3/Lib2)")
p3 <- p3 + coord_cartesian(ylim = c(-0.3,0.3))

p4 <- Plot.Scatter(dat = res.Tg[,c("gene.name", "log2FoldChange", "padj", "gene.length")], 
                   log2FC = log2(1.05), comp.between = "(Lib3/Lib2)")
p4 <- p4 + coord_cartesian(ylim = c(-0.3,0.3))

plots.supp9 <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, res.KO = res.KO, 
                    res.Tg = res.Tg)
save(plots.supp9, file = "../results/SEQC_Brain-Rep_DEGs.RData")

## Num of DEGs Table 
dat <- res.Tg #res.KO
fc <- 1.1
cat("Total Number of DEGs = ", sum(dat$padj < 0.05 & abs(dat$log2FoldChange) > log2(fc), na.rm = TRUE),"\n")
cat("Total Number of Up DEGs = ", sum(dat$padj < 0.05 & dat$log2FoldChange > log2(fc), na.rm = TRUE),"\n")
cat("Total Number of Down DEGs = ", sum(dat$padj < 0.05 & dat$log2FoldChange < -log2(fc), na.rm = TRUE),"\n")

cat("Total Number of Long DEGs = ", sum(dat$padj < 0.05 & abs(dat$log2FoldChange) > log2(fc) &
                                            dat$longGene == "LongGene", na.rm = TRUE),"\n")
cat("Total Number of Long Up DEGs = ", sum(dat$padj < 0.05 & dat$log2FoldChange > log2(fc) &
                                               dat$longGene == "LongGene", na.rm = TRUE),"\n")
cat("Total Number of Long Down DEGs = ", sum(dat$padj < 0.05 & dat$log2FoldChange < -log2(fc) &
                                                 dat$longGene == "LongGene", na.rm = TRUE),"\n")

cat("Total Number of Short DEGs = ", sum(dat$padj < 0.05 & abs(dat$log2FoldChange) > log2(fc) &
                                             dat$longGene == "ShortGene", na.rm = TRUE),"\n")
cat("Total Number of Short Up DEGs = ", sum(dat$padj < 0.05 & dat$log2FoldChange > log2(fc) &
                                                dat$longGene == "ShortGene", na.rm = TRUE),"\n")
cat("Total Number of Short Down DEGs = ", sum(dat$padj < 0.05 & dat$log2FoldChange < -log2(fc) &
                                                  dat$longGene == "ShortGene", na.rm = TRUE),"\n")

