################################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 19th April 2018
#
# Program is used for:
# 1. Nanostring Data analysis
# 2. Comparison of DEGs between Nanostring and Seq (750 genes)
################################################################################

figure5 <- function(){
    ## Run cere_ns() for getting the plots for related to Figure 5 
    load(file = "../results/RTT_cere-SeqvsNS.RData")
    cat("\n\n Fig. 5(A) -- Short Genes\n\n")
    print(cere.seq.ns$figure5[[1]])
    cat("\n\n Fig. 5(B) -- Long Genes\n\n")
    print(cere.seq.ns$figure5[[2]])
    cat("\n\n Fig. 5(C)\n\n")
    print(cere.seq.ns$figure5[[3]])
}

figureS10 <- function(){
    ## Run cere_ns() for getting the plots for related to Figure S10 
    load(file = "../results/RTT_cere-SeqvsNS.RData")
    cat("\n\n Supplementary Fig. 10(A)")
    cat(" -- Mecp2 Cerebellum RNA-seq KO/WT Dataset (Whole Genome)\n\n")
    print(cere.seq.ns$figureS10[[1]]) 
    print(cere.seq.ns$figureS10[[2]]) 
    print(cere.seq.ns$figureS10[[3]]) 
    cat("\n\n Supplementary Fig. 10(B)")
    cat(" -- 750 common genes between RNA-seq and Nanostring\n\n")
    print(cere.seq.ns$figureS10[[4]]) 
    print(cere.seq.ns$figureS10[[5]])
    cat("\n\n Supplementary Fig. 10(C)")
    print(cere.seq.ns$figureS10[[6]]) 
    print(cere.seq.ns$figureS10[[7]])
}

cere_ns <- function() {
    ## Nanostring
    RTT.ns <- read.table("../dat/Mecp2_Seq_NanoString/NS_Mecp2.txt", sep = "\t",
                   stringsAsFactors = FALSE, quote = "", header = TRUE,
                   na.strings = FALSE)
    genotypes <- factor(rep(c("KO", "WT"), each = 3), levels = c("KO", "WT"))
    RTT.ns.norm <- NanoStringNorm(x = RTT.ns, anno = NA, CodeCount = 'sum',
                                  Background = 'mean', 
                                  SampleContent = 'housekeeping.sum',
                                  round.values = FALSE, take.log = FALSE,
                                  return.matrix.of.endogenous.probes = TRUE)
    RTT.ns.norm <- data.frame(RTT.ns.norm)
    p1 <- PCAplot(data = log2(RTT.ns.norm + 1), genotypes = genotypes)
    
    ## DEGs analysis
    dds.KO <- DESeqDataSetFromMatrix( countData = as.matrix(round(RTT.ns.norm)),
                        colData = data.frame(row.names = colnames(RTT.ns.norm),
                                             genotypes = genotypes),
                        design = ~ genotypes)
    sizeFactors(dds.KO) <- rep(1, 6)
    dds.KO <- DESeq(dds.KO, betaPrior = TRUE)
    results(dds.KO, contrast = c("genotypes", "KO", "WT"))
    res.KO <- results(dds.KO, contrast = c("genotypes", "KO", "WT"))
    res.KO <- inner_join(x = data.frame(res.KO) %>% 
                             rownames_to_column(var = "gene.name"),
                         y = data.frame(RTT.ns.norm) %>% 
                             rownames_to_column(var = "gene.name"),
                         by = "gene.name")
    ind.B <- grep(pattern = "KO_", x = colnames(res.KO))
    ind.A <- grep(pattern = "WT_", x = colnames(res.KO))
    res.KO <- logofMeans.between.A.B(dat = res.KO, B.samples = ind.B,
                                     A.samples = ind.A)
    
    ## Cere RNA-Seq annot
    # seq.annot <- read.table("../dat/Mecp2_Seq_NanoString/Cere_Seq_NS_DEGs.txt",
    #                         sep = "\t", stringsAsFactors = FALSE, quote = "",
    #                         header = TRUE)
    seq.plots <- cere_seq()
    seq.annot <- seq.plots$seq.annot
    seq.annot$longGene <- ifelse(seq.annot$gene.length > 100e3, "Long Gene", 
                                 "Short Gene")
    res.KO <- res.KO[, c("gene.name", "log2FoldChange", "logFC.crude", "padj")]
    colnames(res.KO)[2:4] <- paste("NS.", colnames(res.KO)[2:4], sep = "")
    colnames(seq.annot)[2:4] <- paste("RNASeq.", colnames(seq.annot)[2:4], 
                                      sep = "")
    seq.ns <- inner_join(x = res.KO, y = seq.annot, by = "gene.name")
    rownames(seq.ns) <- seq.ns$gene.name
    seq.ns$longGene <- factor(seq.ns$longGene, levels = c("Long Gene",
                                                          "Short Gene"))
    mle.val <- paste("R^2 =", round(summary(
                    lm(formula = RNASeq.logFC.crude ~ NS.logFC.crude,
                       data = seq.ns))$r.squared, 3))
    p2 <- ggplot(data = seq.ns, aes(x=NS.logFC.crude, y=RNASeq.logFC.crude)) +
        geom_point() + geom_abline(color =  "blue") + ylim(c(-1, 1)) + 
        xlim(c(-1, 1)) + annotate("text", x = 0.75, y = -0.75, label = mle.val, 
                                  size = 6, fontface = "bold") + 
        xlab("Nanostring Classical Log2FC") + ylab("RNA-Seq Classical Log2FC") + 
        theme(legend.text = element_text(size = 16, face = "bold"),
              legend.title = element_text(size = 16,colour = "black",
                                          face="bold"), 
              plot.title=element_blank(),
              axis.title = element_text(size = 18,face = "bold"),
              axis.text.x = element_text(size = 16,face = "bold",color="black"),
              axis.text.y = element_text(size = 16,face = "bold",color="black"),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    map.val <- paste("R^2 =", round(summary(
                    lm(formula = RNASeq.log2FoldChange ~ NS.log2FoldChange,
                        data = seq.ns))$r.squared, 3))
    p3 <- ggplot(data = seq.ns, aes(x=NS.log2FoldChange, 
                                    y = RNASeq.log2FoldChange)) + geom_point() + 
        geom_abline(color =  "blue") + ylim(c(-1, 1)) + xlim(c(-1, 1)) +
        annotate("text",x = 0.5,y = -0.75,label = map.val,size = 6,
                 fontface = "bold") + xlab("Nanostring DESeq2 Log2FC") + 
        ylab("RNA-Seq DESeq2 Log2FC") + 
        theme(legend.text = element_text(size = 16, face = "bold"),
              legend.title = element_text(size = 16,colour = "black",
                                          face = "bold"),
            plot.title = element_blank(),
            axis.title = element_text(size = 18, face = "bold"),
            axis.text.x = element_text(size = 16,face = "bold",color = "black"),
            axis.text.y = element_text(size = 16,face = "bold",color = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    seq.ns <- seq.ns[, -c(3,6)]
    seq.ns$logFC.crude <- seq.ns$RNASeq.log2FoldChange
    seq.ns$logFC.crude <- seq.ns$NS.log2FoldChange
    seq.ns$`DEGs Called` <- ifelse( 
                seq.ns$NS.padj < 0.05 & abs(seq.ns$NS.log2FoldChange) > 0 &
                seq.ns$RNASeq.padj < 0.05 & abs(seq.ns$RNASeq.log2FoldChange)>0,
                "Both", ifelse(seq.ns$RNASeq.padj < 0.05 & 
                           abs(seq.ns$RNASeq.log2FoldChange) > 0,"RNA-Seq",
                ifelse(seq.ns$NS.padj < 0.05 &abs(seq.ns$NS.log2FoldChange) > 0,
                        "NanoString","Not Sig")))
    mat <- seq.ns[which(!is.na(seq.ns$`DEGs Called`) & 
                            seq.ns$`DEGs Called` != "Not Sig"), ]
    mat1 <- mat[which(mat$longGene == "Long Gene"), ]
    p5 <- qplot(data = mat1, x = RNASeq.log2FoldChange, y = NS.log2FoldChange,
          col = `DEGs Called`) +
        geom_point(shape = 1.5) + geom_jitter() + 
        geom_abline(colour = "grey", size = 0.5) +
        xlim(c(-1, 1)) + ylim(c(-1, 1)) + ggtitle("Long Genes") +
        geom_hline(aes(yintercept = 0), linetype = "dotted") +
        geom_vline(aes(xintercept = 0), linetype = "dotted") +
        geom_hline(aes(yintercept = log2(1.2)), linetype = "dotted") +
        geom_vline(aes(xintercept = log2(1.2)), linetype = "dotted") +
        geom_hline(aes(yintercept = -log2(1.2)), linetype = "dotted") +
        geom_vline(aes(xintercept = -log2(1.2)), linetype = "dotted") +
        xlab(paste("LogFC RNA-Seq Values")) + ylab(paste("LogFC NanoString Values")) +
        annotate("text", x = 0.85, y = 0.32, label = "log2(1.2)", size = 4,
                 fontface = "bold") +
        annotate("text", x = 0.85, y = -0.32, label = "-log2(1.2)", size = 4,
                 fontface = "bold") +
        theme(legend.text = element_text(size = 14, face = "bold"),
              aspect.ratio = 1, legend.title = element_text(size = 14,
                                                            colour = "black",
                                                            face = "bold"),
            plot.title = element_text(size = 18,colour = "black",face = "bold"),
            axis.title = element_text(size = 18, face = "bold"),
            axis.text.x = element_text(size = 18,face = "bold",color = "black"),
            axis.text.y = element_text(size = 18,face = "bold",color = "black"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    mat2 <- mat[which(mat$longGene == "Short Gene"), ]
    p4 <- qplot(data = mat2, x = RNASeq.log2FoldChange, y = NS.log2FoldChange,
          col = `DEGs Called`) + geom_point(shape = 1) + geom_jitter() + 
        geom_abline(colour = "grey", size = 0.5) + 
        xlim(c(-1, 1)) + ylim(c(-1, 1)) + ggtitle("Short Genes") +
        geom_hline(aes(yintercept = 0), linetype = "dotted") +
        geom_vline(aes(xintercept = 0), linetype = "dotted") +
        geom_hline(aes(yintercept = log2(1.2)), linetype = "dotted") +
        geom_vline(aes(xintercept = log2(1.2)), linetype = "dotted") +
        geom_hline(aes(yintercept = -log2(1.2)), linetype = "dotted") +
        geom_vline(aes(xintercept = -log2(1.2)), linetype = "dotted") +
        xlab(paste("LogFC RNA-Seq Values")) + ylab(paste("LogFC NanoString 
                                                         Values")) +
        annotate("text", x = 0.85, y = 0.32, label = "log2(1.2)",size = 4,
                 fontface = "bold") +
        annotate("text", x = 0.85, y = -0.32, label = "-log2(1.2)", size = 4,
                 fontface = "bold") +
        theme(legend.text = element_text(size = 14, face = "bold"),
              aspect.ratio = 1, legend.title = element_text(size = 14,
                                                            colour = "black",
                                                            face = "bold"),
              plot.title = element_text(size = 18,colour = "black",face="bold"),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(size = 18, face="bold",color="black"),
              axis.text.y = element_text(size = 18,face = "bold",color="black"),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    
    ## Absolute Diff.
    mat <- seq.ns
    mat$abs.diff.log2FC <- apply(mat[, c("RNASeq.log2FoldChange", 
                                         "NS.log2FoldChange")], 1, 
                                 function(r) abs(r[1]) - abs(r[2]))
    up.LongGenes <- sum(mat$abs.diff.log2FC > 0 & mat$gene.length > 100e3, 
                        na.rm = TRUE)
    down.LongGenes <- sum(mat$abs.diff.log2FC < 0 & mat$gene.length > 100e3, 
                          na.rm = TRUE)
    up.ShortGenes <- sum(mat$abs.diff.log2FC > 0 &  mat$gene.length <= 100e3, 
                         na.rm = TRUE)
    down.ShortGenes <- sum(mat$abs.diff.log2FC < 0 & mat$gene.length <= 100e3, 
                           na.rm = TRUE)
    cont.tab <- matrix(data = c(up.LongGenes,down.LongGenes,up.ShortGenes,
                                down.ShortGenes),nrow = 2,ncol = 2)
    rownames(cont.tab) <- c("Up", "Down")
    colnames(cont.tab) <- c("Long.Gene", "Short.Gene")
    print(cont.tab)
    chisq.pval <- paste("Chi Sq. P-Value:", format(round(
                        chisq.test(cont.tab)$p.value, digits = 5),
                        scientific = TRUE))
    mat$longGene <- relevel(x = factor(mat$longGene), ref = "Long Gene")
    p6 <- qplot(data = mat[!is.na(mat$abs.diff.log2FC), ],y = abs.diff.log2FC,
          x = gene.length/1000, xlab = "Gene Length in KB",
          ylab = "Absolute (RNA-Seq - Nanostring) Log2FC") + 
        geom_hline(yintercept = 0) + 
        scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
        theme_grey() + geom_point(aes(colour = longGene)) + ylim(c(-0.4, 0.4)) +
        annotate("text",x = 500,y = 0.25,label = cont.tab[1, 1],size = 7,
                 fontface = "bold") +
        annotate("text",x = 500,y = -0.25,label = cont.tab[2, 1],size = 7,
                 fontface = "bold") +
        annotate("text",x = 1,y = 0.25,label = cont.tab[1, 2],size = 7,
                 fontface = "bold") +
        annotate("text",x = 1,y = -0.25,label = cont.tab[2, 2],size = 7,
                 fontface = "bold") +
        annotate("text",x = 25,y = -0.4,label = chisq.pval,size = 6,
                 fontface = "bold") +
        theme(legend.position = "none",
              plot.title = element_text(size = 14, face = "bold"),
              axis.title = element_text(size = 18, face = "bold"),
              legend.text = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(size=16,face="bold",color = "black"),
              axis.text.y = element_text(size = 18,face="bold",color="black"),
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    cere.seq.ns <- list(figure5 = list(plotA = p4, plotB = p5, plotC = p6),
                        figureS10 = list(plotA1 = seq.plots$p1, 
                                         plotA2 = seq.plots$p2,
                                         plotA3 = seq.plots$p3,
                                         plotB1 = seq.plots$p4, plotB2 = p1, 
                                         plotC1 = p2, plotC2 = p3)) 
    #save(cere.seq.ns, file = "../results/RTT_cere-SeqvsNS.RData")
    return(cere.seq.ns)
}