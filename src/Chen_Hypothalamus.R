###############################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 14th April 2016
#
# Program is used for:
# 1. Long and Short Gene analysis for Lin Chen's Hypothalamus Dataset 
###############################################################################

chen.hypothalamus.dataset <- function(){
  ## dataset
  chen.dataset <- read.table(
      "../dat/DEGS_PreviousStudies/2015-PNAS_Lin/count.table_Lin.txt",
      sep = "\t", stringsAsFactors=F, header = TRUE)
  mm9.ucsc.annot <- read.table(
      file = "../dat-infofiles/annot.files/mm9_ucsc.annot.txt", header = FALSE, 
      sep = "\t",stringsAsFactors = FALSE, quote = "", na.strings = FALSE)
  mm9.ucsc.annot <- mm9.ucsc.annot[,c(1:5)]
  colnames(mm9.ucsc.annot) <- c("ucsc_id", "chr", "strand", "tx_start", 
                                "tx_end")
  mm9.ucsc.annot$gene.length <- mm9.ucsc.annot$tx_end - 
                                    mm9.ucsc.annot$tx_start + 1
  chen.dataset <- inner_join(x = chen.dataset, y = mm9.ucsc.annot, 
                             by = "ucsc_id")
  chen.dataset.exp <- aggregate(.~ gene_name, data = chen.dataset[,c(2:14)], 
                                mean)
  chen.dataset.gl <- aggregate(.~ gene_name, data = chen.dataset[,c("gene_name",
                                                        "gene.length")], max)
  chen.dataset.exp <- inner_join(x = chen.dataset.exp, y = chen.dataset.gl, 
                                 by = "gene_name")
  genotype <- factor(c("WT_C57BL", "KO_C57BL", "WT_C57BL", "KO_C57BL",
                "WT_C57BL", "KO_C57BL", rep(c("Tg_FVB","WT_FVB"), each = 3)), 
                levels = c("KO_C57BL","WT_C57BL","WT_FVB","Tg_FVB"))
  colnames(chen.dataset.exp) <- gsub(pattern = "_count", replacement = "", 
                                     x = colnames(chen.dataset.exp))
  #boxPlot(data = log2(chen.dataset.exp[,c(2:13)]+1), 
  #                 title = "GSE66871/GPL13112", samples = genotype)
  mecp2.id = c("Mecp2","MeCP2")
  rownames(chen.dataset.exp) <- chen.dataset.exp$gene_name
  #GeneExpLevels(data = log2(chen.dataset.exp[,c(2:13)]+1), 
  #             gene.id = mecp2.id, genotypes = genotype, samples = genotype)
  A.samples <- chen.dataset.exp %>% select(contains("wt4ko")) %>% colnames()
  B.samples <- chen.dataset.exp %>% select(contains("ko4ko")) %>% colnames()
  chen.dataset.exp <- logofMeans.between.A.B(dat = chen.dataset.exp, A.samples, 
                                             B.samples)
  chen.dataset.exp$logFC.WT <- apply(X = chen.dataset.exp[,A.samples], 1, 
                                     function(r) (log2((r[3]+1)/(r[1]+1))))
  p1 <- overlay.gabels.plot(mat = chen.dataset.exp[,c("logFC.WT", "logFC.crude",
                                                      "gene.length")], 
                            comp.between1 = "(WT/WT)",comp.between2 = "(KO/WT)")
  cat("\n\n Printing Supp. Figure 2(K) --  Hypothalamus (KO -- RNA-Seq) \n\n")
  print(plot_grid(p1$plot1 + coord_cartesian(ylim = c(-0.5,0.5)), p1$plot2 + 
                      coord_cartesian(ylim = c(0,20)), ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
  A.samples <- chen.dataset.exp %>% select(contains("wt4tg")) %>% colnames()
  B.samples <- chen.dataset.exp %>% select(contains("tg4tg")) %>% colnames()
  chen.dataset.exp <- logofMeans.between.A.B(dat = chen.dataset.exp, A.samples, 
                                             B.samples)
  chen.dataset.exp$logFC.WT <- apply(X = chen.dataset.exp[,A.samples], 1, 
                                     function(r) (log2((r[1]+1)/(r[2]+1))))
  p2 <- overlay.gabels.plot(mat = chen.dataset.exp[,c("logFC.WT", "logFC.crude",
                                                      "gene.length")], 
                            comp.between1 = "(WT/WT)",comp.between2 = "(KO/WT)")
  cat("\n\n Printing Supp. Figure 2(L) -- Hypothalamus (Tg -- RNA-Seq)\n\n")
  print(plot_grid(p2$plot1 + coord_cartesian(ylim = c(-0.25,0.25)), p2$plot2 + 
                      coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
}