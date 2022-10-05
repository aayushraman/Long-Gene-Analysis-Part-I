#############################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 19th April 2018
#
# Program is used for:
# 1. Zylka data analysis
#############################################################################

zylka.analysis <- function(){
  
  ## files
  filePath <- "../dat/Zylka_Dataset/STAR_Analysis/STAR_Counts/"
  file_list <- list.files(path = filePath, pattern = "(.*)_ReadsPerGene.out.tab")
  counts.table <- combiningFilesinDirectory(filePath = filePath, file_list = file_list, n = 2)
  counts.table <- counts.table[-c(1:4),]
  geno.names <- c("Vehicle1", "Vehicle2", "Vehicle3", "Vehicle4", "Vehicle5",
                  "Topotecan1", "Topotecan2", "Topotecan3", "Topotecan4", 
                  "Topotecan5")
  treatments <- factor(c(rep("Untreated_Vehicle", 5), rep("Treated_Topotecan",5)), 
                       levels = c("Untreated_Vehicle", "Treated_Topotecan"))
  mice.strain <- c(rep("Cort.Neurons", 10))
  mm9.annot <- read.table(file = "../dat-infofiles/Id-Symbol.mm9.gencodeIds.txt", 
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                          quote = "", na.strings = FALSE)
  mm9.annot$gene.length <- mm9.annot$end - mm9.annot$start + 1
  counts <- as.matrix(counts.table, colnames = FALSE, row.names =  FALSE)
  colData <- data.frame(treatments, mice.strain, row.names = 
                            colnames(counts.table))
  
  ## DESeq2
  dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ treatments)
  dds <- estimateSizeFactors(dds)
  dat.dds <- counts(dds, normalized = TRUE)
  dat.dds <- data.frame(dat.dds) %>% rownames_to_column()
  colnames(dat.dds)[1] <- "gene.name" 
  dat.dds <- inner_join(x = dat.dds, y = mm9.annot, by = c("gene.name" = "gene"))
  #save(file = "../results/King_RNAseq.RData", dat.dds)
  
  ## Long gene trend
  log2FC.length.WT <- log2FCwithingenotypes(dat = dat.dds[,c(1,2:5,18)])
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.dds, A.samples = c(2:6), 
                                             B.samples = c(7:11))
  log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat3")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude", 
                                                      "gene.length")], 
                              by = "gene.name")
  p1 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], comp.between1 = "(V/V)", 
                            comp.between2 = "(D/V)")
  cat("\n\n Printing Supp. Figure 4 (A) -- Cortical Neurons; 300nM Topotecan treatment (RNA-Seq)\n\n") 
  print(plot_grid(p1$plot1 + coord_cartesian(ylim = c(-2,2)), p1$plot2 + 
                      coord_cartesian(ylim = c(0,40)), ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
  rm(log2FC.length, log2FC.length.WT, log2FC.length.KO)
  library("pd.mouse430.2")
  expt <- "GSE43900_RAW"
  group <- c("G0","G0","G0","G1","G1","G1")
  dat.micro <- read.geo.function(expt = expt, group = group)
  group <- dat.micro$sort.group
  treatments.micro <- factor(c(rep("Topotecan", 3), rep("Vehicle",3)),
                             levels = c("Vehicle", "Topotecan"))
  mice.strain.micro <- c(rep("Cort. Neurons", 6))
  dat1 <- 2^dat.micro$norm.data
  load(file = "../dat-infofiles/mm9.GPL1261.RData")
  colnames(dat1) <- c("GSM1","GSM2","GSM3","GSM4","GSM5","GSM6")
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
  #save(file = "../results/King_Array.RData", dat.annot)
  dat.annot$comp.mat3 <- apply(X = dat.annot[,c(1,5:8)], 1,
        function(r){log2((mean(as.numeric(r[c(2)]))+1)/(as.numeric(r[4])+1))})
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, 
                                             A.samples = c(5:7), 
                                             B.samples = c(2:4))
  log2FC.length <- inner_join(x = dat.annot[,c("gene.name","comp.mat3")],
            y = log2FC.length.KO[,c("gene.name","logFC.crude", "gene.length")], 
            by = "gene.name")
  p1 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], 
                            comp.between1 = "(V/V)", comp.between2 = "(D/V)")
  print(plot_grid(p1$plot1 + coord_cartesian(ylim = c(-2,2)), p1$plot2 + 
                      coord_cartesian(ylim = c(0,40)), ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
  rm(log2FC.length, log2FC.length.KO)
}