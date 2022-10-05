#############################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 18th April 2016
#
# Program is used for:
# 1. Analysis of DEGs genes in Mecp2 studies -- Figure 3
#############################################################################

figure3 <- function(){
    # DEGs from Topotecan Dataset
    load("../results/Topo_DEGS_pub.RData")
    cat("\n\n Fig. 3(A) -- 300nM Topotecan Treatment RNA-seq Datasets \n\n")
    cat("\n\n Cortical Neurons (King et al.) \n\n")
    print(topo.plots$plot1)
    cat("\n\n Cortical Neurons (Mabb et al.) \n\n")
    print(topo.plots$plot2)

    # DEGs from RTT related publication
    #degs.pub <- degs.analysis()
    load("../results/RTT.Degs.RData")
    
    ## Mecp2 related Array (KO/WT)
    cat("\n\n Fig. 3(B) -- Mecp2 related Array (KO/WT) \n\n")
    cat("\n\n Hypothalamus (KO/WT) \n\n")
    print(degs.pub$Maria.KO)
    cat("\n\n Cerebellum (KO/WT) \n\n")
    print(degs.pub$Ben.KO)

    ## Mecp2 related Array (Tg/WT)
    cat("\n\n Fig. 3(C) -- Mecp2 related Array (Tg/WT) \n\n")
    cat("\n\n Hypothalamus (Tg/WT) \n\n")
    print(degs.pub$Maria.Tg)
    cat("\n\n Cerebellum (Tg/WT) \n\n")
    print(degs.pub$Ben.Tg)

    ## Mecp2 related RNA-seq Datasets
    cat("\n\n Fig. 3(D) -- Mecp2 related RNA-seq Datasets \n\n")
    cat("\n\n Hypothalamus (KO/WT) \n\n")
    print(degs.pub$Lin.KO)
    cat("\n\n Hypothalamus (Tg/WT) \n\n")
    print(degs.pub$Lin.Tg)
}

figureS5 <- function(){
    load("../results/RTT.Degs.RData")
    cat("\n\n Supplementary Fig. 5 \n\n")
    cat("\n\n Hippocampus 4 weeks (KO/WT) \n\n")
    print(degs.pub$Baker.KO.4weeks)
    cat("\n\n Hippocampus 9 weeks (KO/WT) \n\n")
    print(degs.pub$Baker.KO.9weeks)
}

mecp2.degs.analysis <- function(){

    # Maria Dataset (Hypothalamus)
    maria.degs.data <- read.table(file = 
    "../dat/DEGS_PreviousStudies/2008-Science_Maria/DegsList.Hypothalamus.txt", 
                                  sep = "\t", header = TRUE, quote = "", 
                                  stringsAsFactors = FALSE)
    maria.degs.data <- maria.degs.data[,-c(10:15)]
    
    ## annot from mgi website
    mm10 <- read.table(file = "../dat-infofiles/annot.files/mm10_hypo.txt", 
                       sep = "\t", header = TRUE, quote = "", 
                       stringsAsFactors = FALSE, fill = TRUE)
    colnames(mm10)[1] <- c("Gene_Symbol")
    mm10$gene.length <- mm10$End - mm10$Start + 1
    mm10.gl <- aggregate(.~ Gene_Symbol, data = mm10[,c("Gene_Symbol",
                                                        "gene.length")], max)
  
    ## 190 genes not present in the analysis (however, most of them 
    ## are loc genes)
    maria.degs.data$Gene_Symbol[!(unique(maria.degs.data$Gene_Symbol) %in% 
                                      unique(mm10.gl$Gene_Symbol))]
    maria.degs.data <- inner_join(x = maria.degs.data, y = mm10.gl,
                                  by = "Gene_Symbol")
    dat.KO <- maria.degs.data[,c(1,6,9,10)]
    maria.KO <- Plot.Scatter(dat = dat.KO, log2FC = 0.2, 
                             comp.between = "(KO/WT)")
    dat.Tg <- maria.degs.data[,c(1,5,8,10)]
    maria.Tg <- Plot.Scatter(dat = dat.Tg, log2FC = 0.2, 
                             comp.between = "(Tg/WT)")
    rm(maria.degs.data,dat.KO,dat.Tg,mm10, mm10.gl)
  
    # Ben-Shachar Dataset (Cerebellum)
    ben.degs.KO.data <- read.table(file = 
  "../dat/DEGS_PreviousStudies/2009-HMG_Ben-Shachar/DegsList.KO.Cerebellum.txt", 
                                     sep = "\t", header = TRUE, quote = "", 
                                     stringsAsFactors = FALSE)
    mm10 <- read.table(file = "../dat-infofiles/annot.files/mm10_cere_KO.txt", 
                       sep = "\t", header = TRUE, quote = "", 
                       stringsAsFactors = FALSE, fill = TRUE)
    colnames(mm10)[1] <- c("Gene_Symbol")
    mm10$gene.length <- mm10$End - mm10$Start + 1
    mm10.gl <- aggregate(.~ Gene_Symbol, data = mm10[,c("Gene_Symbol",
                                                        "gene.length")], max)
  
    ## 85 genes not present in the analysis (however, most of them are loc genes)
    ben.degs.KO.data$Gene_Symbol[!(unique(ben.degs.KO.data$Gene_Symbol) %in% 
                                       unique(mm10.gl$Gene_Symbol))]
    ben.degs.KO.data <- inner_join(x = ben.degs.KO.data, y = mm10.gl, 
                                   by = "Gene_Symbol")
    sum(duplicated(ben.degs.KO.data$Gene_Symbol))
    ind <- which(!duplicated(ben.degs.KO.data$Gene_Symbol))
    ben.degs.KO.data <- ben.degs.KO.data[ind,]
    ben.KO <- Plot.Scatter(dat = ben.degs.KO.data, log2FC = 0.2, 
                            comp.between = "(KO/WT)")
    rm(ben.degs.KO.data)
  
    ## Tg dat
    ben.degs.Tg.data <- read.table(file =
  "../dat/DEGS_PreviousStudies/2009-HMG_Ben-Shachar/DegsList.Tg.Cerebellum.txt",
                                    sep = "\t", header = TRUE, quote = "", 
                                    stringsAsFactors = FALSE)
    mm10 <- read.table(file = "../dat-infofiles/annot.files/mm10_cere_Tg.txt", 
                       sep = "\t", header = TRUE, quote = "", 
                       stringsAsFactors = FALSE, fill = TRUE)
    colnames(mm10)[1] <- c("Gene_Symbol")
    mm10$gene.length <- mm10$End - mm10$Start + 1
    mm10.gl <- aggregate(.~ Gene_Symbol, data = mm10[,c("Gene_Symbol",
                                                        "gene.length")], max)
    ben.degs.Tg.data$Gene_Symbol[!(unique(ben.degs.Tg.data$Gene_Symbol) %in% 
                                       unique(mm10.gl$Gene_Symbol))]
    ben.degs.Tg.data <- inner_join(x = ben.degs.Tg.data, y = mm10.gl, 
                                   by = "Gene_Symbol")
    sum(duplicated(ben.degs.Tg.data$Gene_Symbol))
    ind <- which(!duplicated(ben.degs.Tg.data$Gene_Symbol))
    ben.degs.Tg.data <- ben.degs.Tg.data[ind,]
    ben.Tg <- Plot.Scatter(dat = ben.degs.Tg.data, log2FC = 0.2, 
                           comp.between = "(Tg/WT)")
    rm(ben.degs.Tg.data, mm10.gl, mm10)
  
    # Baker Dataset (Hippocampus)
    dat.KO.4weeks <- read.table(file = 
            "../dat/DEGS_PreviousStudies/2013-Cell_Baker/KO.Hippo4weeks.txt", 
            sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
    mm10 <- read.table(file = "../dat-infofiles/annot.files/mm10_hippo.txt",
                       sep = "\t", header = TRUE, quote = "", 
                       stringsAsFactors = FALSE, fill = TRUE)
    colnames(mm10)[1] <- c("Gene_Symbol")
    mm10$gene.length <- mm10$End - mm10$Start + 1
    mm10.gl <- aggregate(.~ Gene_Symbol, data = mm10[,c("Gene_Symbol",
                                                        "gene.length")], max)
    dat.KO.4weeks$Gene_Symbol[!(unique(dat.KO.4weeks$Gene_Symbol) %in% 
                                    unique(mm10.gl$Gene_Symbol))]
    dat.KO.4weeks <- inner_join(x=dat.KO.4weeks, y=mm10.gl, by="Gene_Symbol")

    ## logFC was calculated as WT/KO in original paper
    ## logFC > 0.1 was used as a cut-off for statistically significance
    dat.KO.4weeks$log2FoldChange <- -dat.KO.4weeks$log2FoldChange
    Baker.KO.4weeks <- Plot.Scatter(dat = dat.KO.4weeks[,c(1,2,4,5)], 
                                    log2FC = 0.1, comp.between = "(KO/WT)")
    rm(dat.KO.4weeks)
  
    ## KO dat -- 9 weeks
    dat.KO.9weeks <- read.table(file = 
            "../dat/DEGS_PreviousStudies/2013-Cell_Baker/KO.Hippo9weeks.txt", 
                                sep = "\t", header = TRUE, quote = "", 
                                stringsAsFactors = FALSE)
    dat.KO.9weeks <- inner_join(x = dat.KO.9weeks, y = mm10.gl, 
                                by = "Gene_Symbol")

    ## logFC was calculated as WT/KO in original paper
    ## logFC > 0.1 was used as a cut-off for statistically significance
    dat.KO.9weeks$log2FoldChange <- -dat.KO.9weeks$log2FoldChange
    Baker.KO.9weeks <- Plot.Scatter(dat = dat.KO.9weeks[,c(1,2,4,5)], 
                                    log2FC = 0.1, comp.between = "(KO/WT)")
    rm(dat.KO.9weeks)
  
    # Lin Dataset (Hypothalamus)
    Lin.KO.WT <- read.table(file =
                    "../dat/DEGS_PreviousStudies/2015-PNAS_Lin/DEGS_WT-KO.txt", 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                    quote = "", na.strings = FALSE)
    mm9.ucsc.annot <- read.table(file =
                            "../dat-infofiles/annot.files/mm9_ucsc.annot.txt", 
                                 header = FALSE, sep = "\t", na.strings = FALSE,
                                 stringsAsFactors = FALSE, quote = "")
    mm9.ucsc.annot <- mm9.ucsc.annot[,c(1:5)]
    colnames(mm9.ucsc.annot) <- c("tx.name", "chr", "strand", "tx_start",
                                  "tx_end")
    mm9.ucsc.annot$gene.length <- mm9.ucsc.annot$tx_end - 
                                  mm9.ucsc.annot$tx_start + 1
    Lin.KO.WT <- inner_join(x = Lin.KO.WT[,c(1,3,5,18)], 
                            y = mm9.ucsc.annot, by = c("gene" = "tx.name"))
    Lin.KO.WT.logFC <- aggregate(.~ name, data = Lin.KO.WT[,c("name", "logFC")],
                                 mean)
    Lin.KO.WT.pval <- aggregate(.~ name, data = Lin.KO.WT[,c("name", "FDR")],
                                min)
    Lin.KO.WT.gl <- aggregate(.~ name, data = Lin.KO.WT[,c("name", 
                                                           "gene.length")], max)
    Lin.KO.WT.logFC <- inner_join(x = Lin.KO.WT.logFC, y = Lin.KO.WT.pval, 
                                  by = "name")
    Lin.KO.WT.logFC <- inner_join(x = Lin.KO.WT.logFC, y = Lin.KO.WT.gl, 
                                  by = "name")
    plot1 <- Plot.Scatter(dat = Lin.KO.WT.logFC, log2FC = 0, 
                          comp.between = "(KO/WT)", pval = 1e-5)
    Lin.KO <- plot1 + coord_cartesian(ylim = c(-3.5,3.25))
    Lin.Tg.WT <- read.table(file = 
                    "../dat/DEGS_PreviousStudies/2015-PNAS_Lin/DEGS_Tg-WT.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                    quote = "", na.strings = FALSE)
    Lin.Tg.WT <- inner_join(x = Lin.Tg.WT[,c(1,3,5,18)], y = mm9.ucsc.annot,
                            by = c("gene" = "tx.name"))
    Lin.Tg.WT.logFC <- aggregate(.~ name, data = Lin.Tg.WT[,c("name", "logFC")], 
                                 mean)
    Lin.Tg.WT.pval <- aggregate(.~ name, data = Lin.Tg.WT[,c("name", "FDR")], 
                                min)
    Lin.Tg.WT.gl <- aggregate(.~ name, data = Lin.Tg.WT[,c("name", 
                                                           "gene.length")], max)
    Lin.Tg.WT.logFC <- inner_join(x = Lin.Tg.WT.logFC, y = Lin.Tg.WT.pval, 
                                  by = "name")
    Lin.Tg.WT.logFC <- inner_join(x = Lin.Tg.WT.logFC, y = Lin.Tg.WT.gl, 
                                  by = "name")
    plot2 <- Plot.Scatter(dat = Lin.Tg.WT.logFC, log2FC = 0, 
                          comp.between = "(Tg/WT)", pval = 1e-5)
    Lin.Tg <- plot2 + coord_cartesian(ylim = c(-2,3.25))
    degs.pub <- list(Maria.KO = maria.KO, Ben.KO = ben.KO, Maria.Tg = maria.Tg,
                     Ben.Tg = ben.Tg, Lin.KO = Lin.KO, Lin.Tg = Lin.Tg,
                     Baker.KO.4weeks = Baker.KO.4weeks, 
                     Baker.KO.9weeks = Baker.KO.9weeks)
#    save(file = "../results/RTT.Degs.RData", degs.pub)
}

topo_degs <- function(){
    
    ## DEGs from King et al.
    King.degs <- read.table(
        file = "../dat/DEGS_PreviousStudies/Topotecan/King_degs-fdr_genes.txt", 
        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    Nabb.degs <- read.table(
        file = "../dat/DEGS_PreviousStudies/Topotecan/Nabb_Degs_Top-Veh.txt", 
        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    gencode <- read.table(
        file = "../dat/DEGS_PreviousStudies/Topotecan/Id-Symbol.mm9.gencodeIds.txt", 
        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    gencode$gene.length <- gencode$end - gencode$start + 1
    King.degs$gene.length <- as.numeric(King.degs$length..kb.)*1000
    King.degs <- King.degs[!is.na(King.degs$gene.length),]
    King.degs$logFC <- King.degs$log2FC
    p1 <- Plot.Scatter(dat = King.degs[,c("gene","logFC", "FDR", "gene.length")], 
                       log2FC = log2(1), comp.between = "(Top/Veh)", pval = 0.05)
    p1 <- p1 + coord_cartesian(ylim = c(-5,5)) 
    
    ## Nabb et al Dataset
    Nabb.degs <- left_join(x = Nabb.degs, y = gencode, by = 
                               c("Genes"="GeneName"))
    Nabb.degs.notpresent.gencode <- read.table(
        file = "../dat/DEGS_PreviousStudies/Topotecan/NabbGenes_notPresent-Gencode.txt", 
        header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    Nabb.degs.notpresent.gencode$gene.length1 <- 
        Nabb.degs.notpresent.gencode$End - Nabb.degs.notpresent.gencode$Start + 1
    Nabb.degs <- left_join(x = Nabb.degs, 
                           y = Nabb.degs.notpresent.gencode[,c("Input", 
                                                               "gene.length1")], 
                           by = c("Genes"="Input"))
    Nabb.degs$gene.length[which(is.na(Nabb.degs$gene.length))] <- 
        Nabb.degs$gene.length1[which(is.na(Nabb.degs$gene.length))]
    p2 <- Plot.Scatter(dat = Nabb.degs[,c("Genes","logFC", "FDR", "gene.length")], 
                       log2FC = log2(1), comp.between = "(Top/Veh)", pval = 0.05)
    p2 <- p2 + coord_cartesian(ylim = c(-6,6)) 
    topo.plots <- list(plot1 = p1, plot2 = p2)
    return(topo.plots)
    #save(topo.plots, file = "../results/Topo_DEGS_pub.RData")
}