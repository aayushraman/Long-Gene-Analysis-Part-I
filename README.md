# Long Gene Analysis (Part I)

Code/Datasets used in the **Apparent bias towards long gene misregulation in MeCP2 syndromes disappears after controlling for baseline variations** project.

-- Dataset: RNA-seq and NanoString datasets that support the findings of this study have been deposited in the Gene Expression Omnibus (GEO) database with the accession codes as follows: [GSE94073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94073), [GSE105047](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105047), and [GSE107399](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107399). GEO accession codes and details of all the public datasets used in this manuscript are available in Supplementary Data 1 table of the paper.  

-- The entire workspace including scripts/processed datasets to reproduce all the results/plots used in this manuscript is available on [zenodo](https://doi.org/10.5281/zenodo.1226607).

The following repo consists of only R scripts used in the project. Please download the dat/dat-infofiles/results directories from the [zenodo workspace](https://doi.org/10.5281/zenodo.1226607) to reporduce the plots. Please contact me directly at aayushraman09@gmail.com if you have any questions/comments regarding the "overlap plots". We are currently working on the R package.

## Example of running overlap plot
```R

## libraries
library(cowplot)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(scales)

## source files from src folder
source("src/overlay.gabel.function.R")
source("src/SEQC_logofMeans.R")
source("src/readGEO.R")

## dat file
counts.table <- read.table(file = "Cere.counts.table_mm10.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(counts.table) <- counts.table$ENSEMBL.ID
counts.table$gene.length <- counts.table$end - counts.table$start + 1

## coldData and dataset with all the gene Ids
genotypes <- factor(c(rep("KO", 3), rep("WT", 3)), levels = c("KO", "WT"))
colData <- data.frame(row.names = colnames(counts.table[, c(2:7)]), genotypes = genotypes)
head(colData)
dim(colData)

## DESeq2
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts.table[, c(2:7)]), colData = colData, design = ~ genotypes)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)

## results from DESeq2
dds <- DESeq(dds, betaPrior = TRUE)
res.KO <- results(dds, contrast = c("genotypes", "KO", "WT"))
res.KO$norm.counts <- counts(dds, normalized = TRUE)
res.KO.annot <- inner_join(x = data.frame(res.KO) %>% tibble::rownames_to_column(var = "ENSEMBL.ID"), y = counts.table[, c(1, 8:14)], by = "ENSEMBL.ID")

## Overlay plots 
res.KO.annot <- logofMeans.between.A.B(res.KO.annot, A.samples = c(11:13), B.samples = c(8:10))
res.KO.annot$comp.mat1 <- apply(X = data.frame(res.KO.annot[,c(11:13)]), 1, function(r) {log2((r[3] + 1)/(r[1] + 1))})
p3 <- overlay.gabels.plot(mat = res.KO.annot[,c("comp.mat1", "logFC.crude", "gene.length")], comp.between1 = "(WT/WT)",comp.between2 = "(KO/WT)")
p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.25, 0.25)),
                p3$plot2 + coord_cartesian(ylim = c(0, 20)), ncol = 1,
                align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
p3
```


## References:

1. AT Raman*, AE Pohodich*, YW Wan, HK Yalamanchili, HY Zoghbi, Z Liu. Apparent bias towards long gene misregulation in MeCP2 syndromes disappears after controlling for baseline variations. [*Nature Communications* (2018)](https://www.nature.com/articles/s41467-018-05627-1) (PMID: 30104565)

2. AT Raman. A research parasite’s perspective on establishing a baseline to avoid errors in secondary analyses. [*GigaScience* (2021)](https://academic.oup.com/gigascience/article/10/3/giab015/6168809) (PMID: 33710326)
