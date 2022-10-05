read.geo.function <- function(expt, group){
  
  ## expt path
  #cat ("Reading expt ", expt,"\n")
  path <- paste("../dat/GEO/", expt, "/", sep="")
  cels <- list.files(path = path, all.files = TRUE, pattern = "CEL")
  #cat("List of cell files: ",cels, "\n")
  
  ## reading the files and normalization
  celfiles <- read.celfiles(filenames = paste(path,cels,sep = ""))
  #cat("Normalizing cell files: ",cels, "\n")
  rma.data <- oligo::rma(celfiles) #, target="core")
  norm.data <- exprs(rma.data)
  
  ## re ordering the samples based on their genotype
  #cat("Reordering the genotypes together: \n")
  norm.data <- norm.data[ , order(group)]
  group <- group[order(group)]
  colnames(norm.data) = gsub(pattern = ".CEL|_(.*).CEL|_(.*).CEL.gz",
                             replacement = "", x = colnames(norm.data))
  #cat("Dimension of the normalized data =",dim(norm.data),"\n")
  
  ## returning the normalized data table and sorted group
  data = list(norm.data = norm.data, sort.group = group)
  return(data)
}

## box plot to check if the data was normalized or not
boxPlot <- function(data, title, samples){
  par(mar=c(2+round(max(nchar(colnames(data)))/2),4,2,1),font=2)
  #pal <- c(rep("#9E0142", each = 3),rep("#5E4FA2", each = 3), rep("#F46D43", each = 4), 
           #rep("#FDAE61", each = 4), rep(brewer.pal(6,"Spectral"), each = 3))
  boxplot(data, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=samples)
  legend("topleft", levels(samples), fill = palette(), bty="n", cex = 0.8, xpd = TRUE)
}

## Gene Expression Level in Microarray data
GeneExpLevels <- function(data, gene.id, genotypes, samples){
  
  ## MeCP2 expression level
  ind = which(rownames(data) == gene.id[1])
  matGene = as.vector(as.matrix(data[ind,]))
  plotDat = data.frame(sampleName = colnames(data),
                       normCounts = matGene,
                       genotype = genotypes,
                       Genotype_Condition = samples)
  print(ggplot(plotDat, aes(x = sampleName, y = normCounts, fill = genotype)) + 
            geom_bar(stat="identity") + ylab("Log2 Normalized Data") + 
            ggtitle(paste("Barplot for Gene ",gene.id[2],sep="")) + theme_grey() +
            facet_grid(. ~ Genotype_Condition,  space = "free", scale="free") + 
            theme(plot.title = element_text(size = 24, face = "bold"),
                  axis.title.y= element_text(size = 20, colour = "black", 
                                             face = "bold"),
                  axis.text.y = element_text(size = 20, colour = "black", 
                                             face = "bold"),
                  axis.text.x = element_blank(), legend.position="none", 
                  axis.ticks.x = element_blank(), axis.title.x =  element_blank(),
                  strip.text = element_text(size = 18, colour = "red", 
                                            face = "bold"))) 
}

## MDS Plot
MDSplot <- function(data, genotypes, conditions){
  mdsDist = cmdscale(d = dist(t(data)), eig = TRUE, k = 2)
  mdsDist = data.frame(genotypes, x = mdsDist$points[,1]/1e4, y = mdsDist$points[,2]/1e4)
  
  ## plot
  if(missing(conditions)){
    print(ggplot(mdsDist, aes(x = x, y = y, color = genotypes)) + 
              geom_point(size = 1) + ## shape = genotypes,
              ylab("MDS Coordinate 2 (x 1e4)") + xlab("MDS Coordinate 1 (x 1e4)") + 
              theme_grey() + theme(legend.text = element_text(size = 18, 
                                                              face = "bold"),
                  legend.title = element_text(size = 18, colour = "black", 
                                              face = "bold"),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 18, face = "bold", 
                                             color = "black"),
                  axis.text.y = element_text(size = 18, face = "bold", 
                                             color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))  
  }
  else{
    print(ggplot(mdsDist, aes(x = x, y = y, shape = conditions, 
                              color = genotypes)) + geom_point(size = 8) +
              ylab("MDS Coordinate 2 (x 1e4)") + xlab("MDS Coordinate 1 (x 1e4)") + 
              theme_grey() + theme(legend.text = element_text(size = 18, 
                                                              face = "bold"),
                  legend.title = element_text(size = 18, colour = "black", 
                                              face = "bold"),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 18, face = "bold", 
                                             color = "black"),
                  axis.text.y = element_text(size = 18, face = "bold", 
                                             color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))  
  }
}

## PCA plot
PCAplot <- function(data, genotypes, conditions){
  ## Calculating PC components
  pcs = prcomp(t(data), center = TRUE)
  percentVar = round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2) 
  if(missing(conditions)){
    #Lib.Id <- genotypes
    print(ggplot(as.data.frame(pcs$x), aes(PC1,PC2, color = genotypes, 
                                           shape = genotypes), 
                 environment = environment()) +
              xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) + 
              geom_point(size = 8) + theme_grey() +
              theme(legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, colour = "black", 
                                              face = "bold"),
                  plot.title = element_blank(),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 16, face = "bold", 
                                             color = "black"),
                  axis.text.y = element_text(size = 16, face = "bold", 
                                             color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
  else{
    print(ggplot(as.data.frame(pcs$x), aes(PC1,PC2, color = genotypes, 
                                           shape = conditions), 
                 environment = environment()) +
            xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) + 
              geom_point(size = 8) + theme_grey() +
            theme(legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, colour = "black", 
                                              face = "bold"),
                  plot.title = element_blank(),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 16, face = "bold", 
                                             color = "black"),
                  axis.text.y = element_text(size = 16, face = "bold", 
                                             color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
}
makeLab = function(x,pc) {
  paste0("PC",pc,": ",x,"% variance")
}

## Scatter Plot for DEGs
Plot.Scatter <- function(dat, log2FC, comp.between, pval = 0.05){
  colnames(dat) <- c("gene.name", "logFC", "adj.P.Val", "gene.length")
  gene.type <- ifelse((dat$adj.P.Val < pval & abs(dat$logFC) > log2FC), 
                      ifelse(dat$gene.length > 100e3, "Long Genes", 
                             "Short Genes"),"Not Stat. Signif.")
  ind <- which(gene.type != "Not Stat. Signif.")
  dat <- dat[ind,]
  gene.type <- gene.type[ind]
  
  ## Contingency Table
  up.LongGenes <- sum(dat$logFC > log2FC & dat$gene.length > 100e3)
  down.LongGenes <- sum(dat$logFC < -log2FC & dat$gene.length > 100e3)
  up.ShortGenes <- sum(dat$logFC > log2FC & dat$gene.length <= 100e3)
  down.ShortGenes <- sum(dat$logFC < -log2FC & dat$gene.length <= 100e3)
  cont.tab <- matrix(data = c(up.LongGenes, down.LongGenes, up.ShortGenes, 
                              down.ShortGenes), nrow = 2, ncol = 2)
  rownames(cont.tab) <- c("Up", "Down")
  colnames(cont.tab) <- c("Long.Gene", "Short.Gene")
  print(cont.tab)
  
  ## qplot
  print(qplot(y = dat$logFC, x = dat$gene.length/1000, colour = gene.type,
              xlab = "Gene Length in KB", ylab = paste("Log2 Fold Change", comp.between)) + 
            scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) + 
            coord_cartesian(ylim = c(-1.5,1.5)) + theme_grey() + 
            annotate("text", x = 500, y=1.5, label= cont.tab[1,1], size=7, fontface="bold") + 
            annotate("text", x = 500, y=-1.5, label= cont.tab[2,1], size=7, fontface="bold") + 
            annotate("text", x = 1, y=1.5, label = cont.tab[1,2], size=7, fontface="bold") +
            annotate("text", x = 1, y=-1.5, label = cont.tab[2,2], size=7, fontface="bold") + 
            theme(plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 18, face = "bold"), legend.position="none", 
                  axis.text.x = element_text(size = 18, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 18, face = "bold", color = "black"),
                  legend.text = element_text(size = 18, face = "bold"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
}

## Scatter Plot with lm line
scatter.lm <- function(dat){  
  r.sq <- paste("R^2 = ",format(summary(lm(gene.length ~ logFC.crude, dat))$r.squared, digits = 2))
  print(summary(lm(gene.length ~ logFC.crude, dat)))
  cat(summary(lm(gene.length ~ logFC.crude, dat))$r.squared,"\n")
  print(qplot(y = logFC.crude, x = gene.length/1000, data = dat) + 
            geom_smooth(method = "lm") + xlab("Gene Length") + 
            ylab("Mean Log2 Fold Change") + 
            scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
            annotate("text", x = 500, y=1.5, label = r.sq, size = 5, 
                     fontface="bold") + 
            theme(plot.title = element_text(size = 18, face = "bold"),
                  axis.title.y = element_text(size = 18, colour = "black", 
                                              face = "bold"),
                    axis.title.x = element_text(size = 18, colour = "black",
                                                face = "bold"),
                    axis.text.y = element_text(size = 18, colour = "black",
                                               face = "bold"),
                    axis.text.x = element_text(size = 18, colour = "black",
                                               face = "bold"),
                    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
}

## box plot with p-value
boxPlot.comp <- function(dat, type, type.dat = "SEQC"){
    if(type.dat == "SEQC"){
        ylab <- paste("Log2 (ß ratio) --",type)
    } else{
        ylab <- paste("Log2",type)    
    }
    pval <- paste("p-value =",format(wilcox.test(logFC.crude ~ longGene, dat)$p.value, 
                                     digits = 2))
    plot <- ggplot(dat, aes(x=longGene, y=logFC.crude)) + theme_classic() + 
            geom_violin(aes(fill = longGene), draw_quantiles = 0.5) +
            geom_jitter(position=position_jitter(width=.1, height=1)) + 
            xlab("") + ylab(ylab) + coord_cartesian(ylim = c(-1,5)) +
            scale_fill_manual(values = c("#00BFC4","#F8766D")) +
            annotate("text", x = 1.5, y = 4, label = pval, size = 5, 
                     fontface="bold") +
            theme(legend.position = "none", 
                  axis.title.y = element_text(size = 18, colour = "black", 
                                              face = "bold"),
                  axis.title.x = element_text(size = 18, colour = "black", 
                                              face = "bold"),
                  axis.text.y = element_text(size = 18, colour = "black", 
                                             face = "bold"),
                  axis.text.x = element_text(size = 18, colour = "black", 
                                             face = "bold"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    p1 <- list(plot = plot, pval = pval)
    return(p1)
}

# ## box plot for Mecp2 dataset
# boxPlot.comp.mecp2 <- function(dat, type){
#     pval <- format(t.test(logFC.crude ~ longGene, dat, 
#                           alternative = "greater")$p.value, digits = 2)
#     pval <- paste("p-value =",pval)
#     plot <- ggplot(dat, aes(x=longGene, y=logFC.crude)) + #geom_boxplot(outlier.shape=NA) +
#               geom_violin(aes(fill = longGene), draw_quantiles = 0.5) +
#               geom_jitter(position=position_jitter(width=.1, height=1)) +
#               xlab("") + ylab(paste("Log2",type)) + theme_grey() +
#               theme(legend.position = "none",axis.title.y = 
#                         element_text(size = 18, colour = "black",face = "bold"),
#                     axis.title.x = element_text(size = 18, colour = "black", 
#                                                 face = "bold"),
#                     axis.text.y = element_text(size = 18, colour = "black", 
#                                                face = "bold"),
#                     axis.text.x = element_text(size = 18, colour = "black", 
#                                                face = "bold"),
#                     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
#     p1 <- list(plot = plot, pval = pval)
#     return(p1)
# }