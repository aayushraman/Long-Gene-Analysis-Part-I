## Deng Dataset GSE6955
deng.human.analysis <- function(){
  cat("Running Deng's Mecp2 Dataset -- Human Dataset\n") 
  Sys.sleep(2)
  
  ## reading the cel files
  expt <- "GSE6955_RAW"
  group <- c("G01","G02","G01","G02","G01","G02")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  #head(dat$norm.data)
  dat$norm.data <- 2^dat$norm.data
  
  ## for plots
  genotypes <- factor(c(rep("WT",3),rep("KO",3)),levels = c("KO", "WT"))
  mice.strain <- c(rep("Human Samples",6))
  sample.info <- factor(c("Normal_2years","Normal_5years","Normal_10years","RTT_2years","RTT_6years","RTT_8years"),
                        levels = c("Normal_2years","RTT_2years","Normal_5years","RTT_6years", "Normal_10years","RTT_8years"))
  ## plots
  # title <- "GSE6955/GPL8300 samples"
  # mecp2.id = c("34355_at","MeCP2")
  # boxPlot(data = log2(dat$norm.data+1), title = title, samples = sample.info)
  # GeneExpLevels(data = log2(dat$norm.data+1), gene.id = mecp2.id, genotypes = genotypes, 
  #               samples = sample.info)
  # PCAplot(data = log2(dat$norm.data+1), genotypes, sample.info)
  
  ## merging anot file without gene length
  data.annot <- read.table(file = "../dat-infofiles/HG_U95Av2.hg38.txt",
                           sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, na.strings = FALSE, fill = TRUE)
  data.annot <- merge(y = data.annot, x = dat$norm.data, by.x = "row.names", by.y = "Probeset_ID")
  gene.length <- data.annot[,c(8:14)]
  gene.length <- gene.length[!duplicated(gene.length),]
  gene.length$gene.length <- gene.length$Gene_End -gene.length$Gene_Start + 1
  
  ## gene length
  dat.agg <- aggregate(. ~ Gene_Name, data = data.annot[,c(2:7,11)], mean)
  dat.agg.length <- merge(x = dat.agg, y = gene.length, by.y = "Gene_Name", by.x = "Gene_Name")
  dat.agg.length <- dat.agg.length[!duplicated(dat.agg.length$Gene_Name),]
  rownames(dat.agg.length) <- dat.agg.length[,1]
  dat.agg.length <- dat.agg.length[,-1]
  
  ## Gabel Plot -- Whole
  WT <- which(group == "G01")
  KO <- which(group == "G02")
  log2FC.length <- logofMeans.between.A.B(dat.agg.length, B.samples = KO, A.samples = WT)
  #gabels.plot(log2FC.length[,c("logFC.crude", "gene.length")])
  
  ## One Sample plot
  dat.agg.length$log2.KO.WT1 <- apply(dat.agg.length[,c(1,4)], 1, function(r){log2((r[c(2)])/(r[c(1)]))})
  #plot1 <- gabels.plot(dat.agg.length[,c("log2.KO.WT1", "gene.length")]) 
  #plot1$plot + coord_cartesian(ylim = c(-0.25,0.15))
  
  dat.agg.length$log2.KO.WT2 <- apply(dat.agg.length[,c(2,5)], 1, function(r){log2((r[c(2)])/(r[c(1)]))})
  #plot2 <- gabels.plot(dat.agg.length[,c("log2.KO.WT2", "gene.length")])
  #plot2$plot + coord_cartesian(ylim = c(-0.15,0.15))
  
  dat.agg.length$log2.KO.WT3 <- apply(dat.agg.length[,c(3,6)], 1, function(r){log2((r[c(2)])/(r[c(1)]))})
  #plot3 <- gabels.plot(dat.agg.length[,c("log2.KO.WT3", "gene.length")])
  #plot3$plot + coord_cartesian(ylim = c(-0.15,0.15))
  
  ## line plots
  dat <- inner_join(x = dat.agg.length %>% add_rownames(), y = log2FC.length %>% add_rownames(), 
                    by = "rowname")
  dat <- dat[,c("logFC.crude", "log2.KO.WT1", "log2.KO.WT2", "log2.KO.WT3", "rowname", "gene.length.x")]
  dat <- dat[order(dat$gene.length.x),]
  dat$gene.length.x <- dat$gene.length.x/1000
  mean.points <- data.frame()
  bin.size <- 200
  shift.size <- 40
  num.bins <- round((dim(dat.agg.length)[1]-bin.size)/shift.size)+1
  
  ## taking the mean of log2FC and genomic length
  for(i in 0:num.bins){
    start <- i*shift.size+1
    end <- start + bin.size-1
    
    ## if the start exceeds total number of genes
    if ((start > dim(dat)[1])) break;
    
    ## if the last bin exceeds the number of genes available
    if(end > dim(dat)[1]){
      end <- dim(dat)[1]
    }
    mat1 <- dat$logFC.crude[start:end]
    mat2 <- dat$log2.KO.WT1[start:end]
    mat3 <- dat$log2.KO.WT2[start:end]
    mat4 <- dat$log2.KO.WT3[start:end]
    mat.length <- mean(dat$gene.length.x[start:end])
    mat.mean1 <- mean(mat1, na.rm = TRUE)
    mat.mean2 <- mean(mat2, na.rm = TRUE)
    mat.mean3 <- mean(mat3, na.rm = TRUE)
    mat.mean4 <- mean(mat4, na.rm = TRUE)
    bin.width <- end-start+1
    
    ## mat means
    #cat("Bin",i+1,": ",start,"-",end,"\t",mat.length,"\t",mat.mean1,"\t",mat.mean2,"\t",mat.mean3,"\t",mat.mean4,"\n",sep="")
    mat.mean <- data.frame(mat.mean1, mat.mean2, mat.mean3, mat.mean4, mat.length)
    mean.points <- rbind(mean.points, mat.mean)
    
    ## end exceeds total number of genes
    if (end == dim(dat)[1]) break;
  }
  
  ## overlay line plot -- 1
  cat("\n\n Printing Figure 4(A)\n\n")
  ind <- mean.points$mat.length >=1 & mean.points$mat.length <=1000
  mean.points <- mean.points[ind, ]
  plot1 <- ggplot(data = mean.points[,c(1,2,5)], aes(x = mat.length)) + 
    geom_line(aes(y = mean.points$mat.mean1, color = colnames(mean.points)[1]), size=1) + 
    geom_line(aes(y = mean.points$mat.mean2, color = colnames(mean.points)[2]), size=1) + 
    #geom_line(aes(y = mean.points$mat.mean3, color = colnames(mean.points)[3]), size=1) +
    #geom_line(aes(y = mean.points$mat.mean4, color = colnames(mean.points)[4]), size=1) + 
    scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
    xlab("Mean Gene Length in KB") + ylab("Mean Log2 Fold Change (RTT/WT)") + theme_grey() +
    scale_color_manual(name = "Age", values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),
                       labels = c("Whole Dataset", "2/4 years", "5 years", "8 years")) +
    theme(plot.title = element_blank(), legend.title = element_text(size = 18, face = "bold"),
          legend.position = c(.85,.2), legend.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 24, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "black"),
          axis.text.y = element_text(size = 22, face = "bold", color = "black"),
          plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))
  print(plot1)
  
  ## overlay line plot -- 2
  cat("\n\n Printing Figure 4(B)\n\n")
  plot2 <- ggplot(data = mean.points[,c(3:5)], aes(x = mat.length)) + 
    geom_line(aes(y = mean.points$mat.mean3, color = colnames(mean.points)[1]), size=1) + 
    geom_line(aes(y = mean.points$mat.mean4, color = colnames(mean.points)[2]), size=1) + 
    scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) + coord_cartesian(ylim = c(-0.15,0.15)) +
    xlab("Mean Gene Length in KB") + ylab("Mean Log2 Fold Change (RTT/Normal)") + theme_grey() +
    scale_color_manual(name = "Age", values = c("#4DAF4A","#984EA3"), labels = c("5 years", "8 years")) +
    theme(plot.title = element_blank(), legend.title = element_text(size = 18, face = "bold"),
          legend.position = c(.85,.2), legend.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 24, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "black"),
          axis.text.y = element_text(size = 22, face = "bold", color = "black"),
          plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))
  print(plot2)
}