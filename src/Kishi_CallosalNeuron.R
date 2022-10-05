## Callosal Neurons RTT Dataset
kishi.callosalneurons.dataanalysis <- function(){
  cat("Running Kishi's Callosal Projection Neurons Dataset\n")
  library("pd.mouse430.2")
  load(file = "../dat-infofiles/mm9.GPL1261.RData")
  
  experiment <- pmid2MIAME("26821816")
  experiment@name <- "Kishi, Noriyuki"
  experiment@lab <- "Jeffrey Macklis Lab, Harvard"
  experiment@contact <- "Jeffrey Macklis Lab or for raw dataset: GSE50225"
  experiment@url <- "https://www.ncbi.nlm.nih.gov/pubmed/26821816"
  annotation <- "GPL1261/Affymetrix Mouse Genome 430 2.0 Array -- NetAffx mm9"
  expt <- "GSE50225_RAW"
  group <- c("G1","G1","G1","G2","G2","G2")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  mice.strain <- c(rep("C57BL",6))
  genotypes <- factor(c(rep("WT", 3), rep("KO",3)),levels = c("KO", "WT"))
  sample.info <- factor(c(rep("WT_C57BL", 3), rep("KO_C57BL", 3)),
                        levels = c("KO_C57BL", "WT_C57BL"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(
                                sample.name = colnames(dat$norm.data),
                                sample.type = sample.info, 
                                mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData, 
                        experimentData = experiment, annotation = annotation)
  box.plot <- boxPlot(data = dat$norm.data, samples = sample.info, title = "")
  mecp2.id <- c("1460246_at","MeCP2")
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, eset$sample.type)
  dat1 <- 2^(exprs(eset))
  mds.plot <- MDSplot(data = dat1, genotypes = genotypes, conditions = mice.strain)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
  RTT.callosal.neurons <- list(dat.annot = dat.annot, box.plot = box.plot, 
                                mds.plot = mds.plot, mecp2.plot = mecp2.plot,
                                eset = eset)
  
  ## Kishi's Dataset
  log2FC.length.WT <- dat.annot[,c(1:4,8)]
  log2FC.length.WT$comp.mat <- apply(dat.annot[,c(2:4)], 1, 
                                     function(r){log2((r[1]+1)/(r[2]+1))})
  log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, 
                                             A.samples = c(2:4), 
                                             B.samples = c(5:7))
  log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude", 
                                                      "gene.length")], 
                              by = "gene.name")
  p3 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], 
                            comp.between1 = "(WT/WT)", 
                            comp.between2 = "(KO/WT)") 
  p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.25,0.25)), p3$plot2 + 
                      coord_cartesian(ylim = c(0,20)),ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  RTT.callosal.neurons <- list(plot = p3)
  #save(file = "../results/RTT_callosal.neurons.RData", RTT.callosal.neurons)
  return(RTT.callosal.neurons)
}