#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 14th April 2018
#
# Program is used for:
# 1. Long and Short Gene analysis -- Supp 4,5 and 6
#############################################################################

sugino.neuron.dataset <- function(){
  
  ## GSE8720_RAW
  library(RColorBrewer)
  cat("Running 	Sugino's Mecp2 Dataset -- First Long gene paper\n") 
  load("../dat-infofiles/mm9.GPL1261.RData")
  
  ## reading the cel files
  expt <- "GSE8720_RAW"
  group <- c(rep("G01",3), rep("G02",3), rep("G03",4), rep("G04",4), rep("G05",3), rep("G06",3), 
             rep("G07",3), rep("G08",3), rep("G09",3), rep("G10",3))
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  
  ## for plots
  genotypes <- factor(c(rep("WT", 3), rep("KO",3), rep("WT", 4), rep("KO",4), rep("WT", 3), 
                        rep("KO",3), rep("WT", 3), rep("KO",3), rep("WT", 3), rep("KO",3)),levels = c("KO", "WT"))
  mice.strain <- c(rep("LC_YT", 6), rep("G42-M1",8), rep("YPFH-M1",6), rep("LC", 6), rep("G42-Cere", 6))
  sample.info <- factor(c(rep("TH-LC_WT_YT", 3), rep("TH-LC_KO_YT", 3), rep("G42-M1_WT", 4), rep("G42-M1_KO", 4), 
                          rep("YPFH-M1_WT", 3), rep("YPFH-M1_KO", 3), rep("TH-LC_WT", 3), rep("TH-LC_KO", 3), 
                          rep("G42-Cere_WT",3), rep("G42-Cere_KO",3)),
                        levels = c("TH-LC_WT_YT","TH-LC_KO_YT","G42-M1_WT","G42-M1_KO","YPFH-M1_WT","YPFH-M1_KO","TH-LC_WT",
                                   "TH-LC_KO","G42-Cere_WT","G42-Cere_KO"))
  ## box plot
  title <- "GSE8720/GPL1261 samples"
  colnames(dat$norm.data) <- gsub(pattern = ".gz", replacement = "", x = colnames(dat$norm.data))
  #boxPlot(data = dat$norm.data, title = title, samples = sample.info)
  
  ## all the plots
  #mecp2.id = c("1438538_at","MeCP2")
  #GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, genotypes = genotypes, samples = sample.info)
  # mecp2.id = c("1438930_s_at","MeCP2")
  # GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, genotypes = genotypes, samples = sample.info)
  # mecp2.id = c("1449266_at","MeCP2")
  # GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, genotypes = genotypes, samples = sample.info)
  # mecp2.id = c("1460246_at","MeCP2")
  # GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, genotypes = genotypes, samples = sample.info)
  #PCAplot(data = dat$norm.data, genotypes, mice.strain)
  #MDSplot(data = dat$norm.data, genotypes, mice.strain)
  
  ## merging 
  dat1 <- 2^(dat$norm.data)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
  
  ## overlay gabel plot for WT/WT and KO/WT LC-Y
  dat.annot$logFC.LC_Y.WT <- apply(dat.annot[,c(2:4)], 1, function(r) {log2(r[2]/r[3])})
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, A.samples = c(2:4), B.samples = c(5:7))
  log2FC.length <- inner_join(x = dat.annot[,c("gene.name","logFC.LC_Y.WT")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude", "gene.length")], by = "gene.name")
  p1 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", comp.between2 = "(KO/WT)")
  cat("\n\n Printing Figure Supp. Figure 2 -- Locus Coeruleus (P22)\n\n")
  print(plot_grid(p1$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p1$plot2 + coord_cartesian(ylim = c(0,10)), 
            ncol = 1, align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
  p1 <- plot_grid(p1$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p1$plot2 + coord_cartesian(ylim = c(0,10)), 
                  ncol = 1, align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  rm(log2FC.length, log2FC.length.KO)
  
  ## overlay gabel plot for WT/WT and KO/WT LC
  dat.annot$logFC.LC.WT <- apply(dat.annot[,c(22:24)], 1, function(r) {log2(r[2]/r[1])})
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, A.samples = c(22:24), B.samples = c(25:27))
  log2FC.length <- inner_join(x = dat.annot[,c("gene.name","logFC.LC.WT")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude", "gene.length")], by = "gene.name")
  p2 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", comp.between2 = "(KO/WT)")
  cat("\n\n Printing Supp. Figure 2 --  Locus Coeruleus\n\n")
  print(plot_grid(p2$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p2$plot2 + coord_cartesian(ylim = c(0,10)), 
            ncol = 1, align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
  p2 <- plot_grid(p2$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p2$plot2 + 
                      coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') +
        theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  rm(log2FC.length, log2FC.length.KO)
  
  ## overlay gabel plot for WT/WT and KO/WT FS
  dat.annot$logFC.FS.WT <- apply(dat.annot[,c(8:11)], 1, function(r) {log2(mean(r[c(2,3)])/mean(r[c(1,4)]))})
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, A.samples = c(8:11), B.samples = c(12:15))
  log2FC.length <- inner_join(x = dat.annot[,c("gene.name","logFC.FS.WT")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude", "gene.length")], by = "gene.name")
  p3 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", comp.between2 = "(KO/WT)")
  cat("\n\n Printing Supp. Figure 2 -- FS, Motor Cortex\n\n")
  p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p3$plot2 + coord_cartesian(ylim = c(0,10)), 
            ncol = 1, align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  rm(log2FC.length, log2FC.length.KO)
  
  ## overlay gabel plot for WT/WT and KO/WT PK
  dat.annot$logFC.PK.WT <- apply(dat.annot[,c(28:30)], 1, function(r) {log2(r[3]/r[2])})
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, A.samples = c(28:30), B.samples = c(31:33))
  log2FC.length <- inner_join(x = dat.annot[,c("gene.name","logFC.PK.WT")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude", "gene.length")], by = "gene.name")
  p4 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", comp.between2 = "(KO/WT)")
  cat("\n\n Printing Supp. Figure 2 -- Purkinje Cells, Cerebellum\n\n")
  print(plot_grid(p4$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p4$plot2 + coord_cartesian(ylim = c(0,10)), 
            ncol = 1, align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
  p4 <- plot_grid(p4$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p4$plot2 + 
                      coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
        theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  rm(log2FC.length, log2FC.length.KO)
  
  ## overlay gabel plot for WT/WT and KO/WT Pyramidal neurons
  dat.annot$logFC.PN.WT <- apply(dat.annot[,c(16:18)], 1, function(r) {log2(r[2]/r[1])})
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, A.samples = c(16:18), B.samples = c(19:21))
  log2FC.length <- inner_join(x = dat.annot[,c("gene.name","logFC.PN.WT")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude", "gene.length")], by = "gene.name")
  p5 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", comp.between2 = "(KO/WT)")
  cat("\n\n Printing Supp. Figure 2 -- PN, Motor Cortex\n\n")
  print(plot_grid(p5$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p5$plot2 + 
                      coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
  p5 <- plot_grid(p5$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p5$plot2 + 
                      coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  rm(log2FC.length, log2FC.length.KO)
  Sugino_Mecp2 <- list(p1,p2,p3,p4,p5)
  save(Sugino_Mecp2, file = "../results/RTT_Sugino_Neurons.RData")
}