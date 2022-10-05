## Maria Hypothalamus
maria.hypothalamus.dataanalysis <- function(){
  cat("Running Maria's Hypothalamus Dataset\n") 
  load("../dat-infofiles/mm9.GPL6193.RData")
  experiment <- new("MIAME", name="RTT Hypothalamus Dataset", 
                    lab="Huda Zoghi Lab, BCM", 
                    title="MeCP2, a Key Contributor to Neurological Disease, 
                    Activates and Represses Transcription")
  annotation <- "GPL6193/Affymetrix Mouse Exon 1.0 ST Array"

  ## reading and normalizing the cel files
  library(pd.moex.1.0.st.v1)
  expt <- "GSE11150_RAW"
  group <- c("G0","G0","G3","G3","G3","G0","G0","G3","G2","G1","G2","G1","G1",
             "G1","G2","G2")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  mice.strain <- c(rep("FVB", 4), rep("C57BL",8), rep("FVB", 4))
  genotypes <- factor(c(rep("WT", 4), rep("WT",4), rep("KO", 4), rep("Tg3",4)),
                      levels = c("KO", "WT", "Tg3"))
  sample.info <- factor(c(rep("WT_FVB", 4), rep("WT_C57BL", 4), 
                          rep("KO_C57BL", 4), rep("Tg3_FVB", 4)),
                        levels = c("KO_C57BL", "WT_C57BL", "WT_FVB", "Tg3_FVB"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(
                                    sample.name = colnames(dat$norm.data),
                                    sample.type = sample.info,
                                    mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset.hypo <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData,
                             experimentData = experiment, annotation = annotation)
  boxPlot(data = dat$norm.data, title = "", samples = sample.info)
  mecp2.id <- c("7017610","MeCP2")
  mecp2.plot <- GeneExpLevels(data = exprs(eset.hypo), gene.id = mecp2.id, 
                              genotypes = genotypes, eset.hypo$sample.type)
  dat1 <- 2^(exprs(eset.hypo))
  mds.plot <- MDSplot(data = dat1, genotypes = genotypes, 
                      conditions = mice.strain)
  
  ## merging the exp mat and annot
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
# RTT.hypo <- list(dat.annot = dat.annot, mds.plot = mds.plot, 
#                  mecp2.plot = mecp2.plot, eset = eset.hypo)
# save(file = "../results/RTT_hypo.RData",RTT.hypo)
  return(list(dat.annot = dat.annot, mds.plot = mds.plot, mecp2.plot = 
                  mecp2.plot, eset = eset.hypo))
}

## Samaco Anygdala
samaco.amygdala.dataanalysis <- function(){
  cat("Running Samaco's Amygdala Dataset\n")
  load("../dat-infofiles/mm9.GPL6193.RData")
  experiment <- new("MIAME", name="RTT Amygdala Dataset", 
                    lab="Huda Zoghi Lab, BCM", 
    title="Crh and Oprm1 mediate anxiety-related behavior and social approach 
    in a mouse model of MECP2 duplication syndrome.")
  annotation <- "GPL6193/Affymetrix Mouse Exon 1.0 ST Array"
  library(pd.moex.1.0.st.v1)
  expt <- "GSE33457_RAW"
  group <- c("G0","G1","G1","G1","G1","G0","G0","G0","G2","G3","G3","G2","G3",
             "G2","G3","G3","G2","G2")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  genotypes <- factor(c(rep("WT", 4), rep("KO",4), rep("Tg", 5), rep("WT",5)),
                      levels = c("KO", "WT", "Tg"))
  mice.strain <- c(rep("B6", 8), rep("FVB", 10))
  sample.info <- factor(c(rep("WT_B6", 4), rep("KO_B6", 4), 
                          rep("Tg3_FVB", 5), rep("WT_FVB", 5)),
                        levels = c("KO_B6", "WT_B6", "WT_FVB", "Tg3_FVB"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(sample.name = 
                                                    colnames(dat$norm.data),
                                                    sample.type = sample.info,
                                                    mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset.amygdala <- ExpressionSet(assayData=dat$norm.data,phenoData = phenoData,
                             experimentData=experiment, annotation=annotation)
  mecp2.id = c("7017610","MeCP2") 
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, sample = sample.info)
  dat1 <- 2^(exprs(eset.amygdala))
  mds.plot <- MDSplot(data = dat1, genotypes = genotypes, 
                      conditions = mice.strain)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
#  RTT.amygdala <- list(dat.annot = dat.annot, mds.plot = mds.plot, 
#                       mecp2.plot = mecp2.plot, eset = eset.amygdala)
#  save(file = "../results/RTT_amygdala.RData",RTT.amygdala)
  return(list(dat.annot = dat.annot, mds.plot = mds.plot, 
              mecp2.plot = mecp2.plot, eset = eset.amygdala))
}  

## Striatum Amygdala
zhao.striatum.dataanalysis <- function(){
  cat("Running Zhao's Striatum Dataset\n")
  load("../dat-infofiles/mm9.GPL6193.RData")
  experiment <- new("MIAME", name="RTT Striatum Dataset", lab="Zhou Lab, UPenn", 
                    title="Loss of MeCP2 function is associated with distinct 
                    gene expression changes in the striatum.")
  annotation <- "GPL6096/Affymetrix Mouse Exon 1.0 ST Array (Gene Version)"
  expt <- "GSE42895_RAW"
  group <- c("G0","G0","G0","G0","G0","G1","G1","G1","G1","G1","G2","G2","G2",
             "G2","G2","G3","G3","G3","G3","G3")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  genotypes <- factor(c(rep("WT", 5), rep("KO",5), rep("WT", 5), rep("KO",5)),
                      levels = c("KO", "WT", "Tg"))
  mice.strain <- c(rep("Striatum", 10), rep("Liver",10))
  sample.info <- factor(c(rep("WT_Striatum", 5), rep("KO_Striatum", 5), 
                          rep("WT_Liver", 5), rep("KO_Liver", 5)),
                        levels = c("KO_Striatum", "WT_Striatum", "KO_Liver",
                                   "WT_Liver"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(sample.name = 
                                                    colnames(dat$norm.data),
                                                    sample.type = sample.info,
                                                    mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset.striatum <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData,
                                 experimentData = experiment, 
                                 annotation = annotation)
  mecp2.id = c("7017610","MeCP2") 
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, sample = sample.info)
  dat1 <- 2^(exprs(eset.striatum))
  mds.plot <- MDSplot(data = dat1, genotypes = genotypes, 
                      conditions = mice.strain)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
# RTT.striatum <- list(dat.annot = dat.annot, mds.plot = mds.plot, 
#                      mecp2.plot = mecp2.plot, eset = eset.striatum)
# save(file = "../results/RTT_striatum.RData",RTT.striatum)
  return(list(dat.annot = dat.annot, mds.plot = mds.plot, mecp2.plot = 
                  mecp2.plot, eset = eset.striatum))
}

## Baker Hippocampus
baker.hippocampus.dataanalysis <- function(){
  cat("Running Baker Hippocampus Dataset\n")
  load("../dat-infofiles/mm9.GPL6246.RData")
  experiment <- new("MIAME", name="RTT Hippocampus Dataset", 
                    lab="Huda Zoghi Lab, BCM", 
                    title="An AT-Hook Domain in MeCP2 Determines the Clinical 
                    Course of Rett Syndrome and Related Disorders.")
  annotation <- "GPL6246/Affymetrix Mouse Gene 1.0 ST Array"
  library(pd.mogene.1.0.st.v1)
  expt <- "GSE42987_RAW"
  group <- c("G0","G0","G1","G1","G2","G2","G3","G3","G3","G3","G2","G2",
             "G1","G1","G0","G0","G4","G4","G5","G5","G6","G6","G7","G7",
             "G7","G7","G6","G6","G5","G5","G4","G4")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  genotypes <- factor(c(rep("WT", 4), rep("KO",4), rep("R270X",4), 
                        rep("G273X",4), rep("WT", 4), rep("KO",4), 
                        rep("R270X",4), rep("G273X",4)),
                      levels = c("KO", "R270X", "G273X", "WT"))
  mice.strain <- c(rep("4 weeks", 16), rep("9 weeks",16))
  sample.info <- factor(c(rep("WT_4weeks", 4), rep("KO_4weeks",4), 
                          rep("R270X_4weeks", 4), rep("G273X_4weeks",4), 
                          rep("WT_9weeks", 4), rep("KO_9weeks",4),
                          rep("R270X_9weeks", 4), rep("G273X_9weeks",4)),
                          levels = c("WT_4weeks", "KO_4weeks", "R270X_4weeks", 
                                     "G273X_4weeks", "WT_9weeks", "KO_9weeks", 
                                     "R270X_9weeks", "G273X_9weeks"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(
                                        sample.name = colnames(dat$norm.data),
                                        sample.type = sample.info,
                                        mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset.hippo <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData,
                                 experimentData = experiment, 
                              annotation = annotation)
  mecp2.id <- c("10605247","MeCP2") 
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, sample = sample.info)
  dat1 <- 2^(exprs(eset.hippo))
  mds.plot <- MDSplot(data = dat1, genotypes = genotypes, 
                      conditions = mice.strain)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
  # RTT.hippo <- list(dat.annot = dat.annot, mds.plot = mds.plot, 
  #                   mecp2.plot = mecp2.plot, eset = eset.hippo)
  # save(file = "../results/RTT_hippo.RData",RTT.hippo)
  return(list(dat.annot = dat.annot, mds.plot = mds.plot, 
              mecp2.plot = mecp2.plot, eset = eset.hippo))
}

## Deng Human RTT Dataset
deng.human.dataanalysis <- function(){
  cat("Running 	Deng's Mecp2 Dataset -- Human Dataset\n") 
  experiment <- new("MIAME", name="Deng Human Dataset", 
                    lab="S.R. Ojeda Lab, OHSU", 
                    title="FXYD1 is an MeCP2 target gene overexpressed in the 
                    brains of Rett syndrome patients and Mecp2-null mice.")
  annotation <- "GPL8300/Affymetrix Human Genome U95 Version 2 Array"
  library(pd.hg.u95av2)
  expt <- "GSE6955_RAW"
  group <- c("G01","G02","G01","G02","G01","G02")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  group
  head(dat$norm.data)
  genotypes <- factor(c(rep("WT",3),rep("KO",3)),levels = c("KO", "WT"))
  mice.strain <- c(rep("Human Samples",6))
  sample.info <- factor(c("Normal_2years","Normal_5years","Normal_10years",
                          "RTT_2years","RTT_6years","RTT_8years"),
                        levels = c("Normal_2years","RTT_2years","Normal_5years",
                                   "RTT_6years", "Normal_10years","RTT_8years"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(
                                          sample.name = colnames(dat$norm.data),
                                          sample.type = sample.info,
                                          mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData,
                              experimentData = experiment, 
                        annotation = annotation)
  mecp2.id <- c("34355_at","MeCP2")
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, samples = sample.info)
  mds.plot <- MDSplot(data = 2^(exprs(eset)), genotypes, sample.info)
  featureData <- read.table(file = "../dat-infofiles/HG_U95Av2.hg38.txt",
                           sep = "\t", header = TRUE, quote = "", 
                           stringsAsFactors = FALSE, na.strings = FALSE, 
                           fill = TRUE)
  data.annot <- inner_join(x = data.frame(2^dat$norm.data) %>% 
                               rownames_to_column(), y = featureData, 
                           by = c("rowname"= "Probeset_ID"))
  gene.length <- data.annot %>% dplyr::select(contains("Gene_"))
  gene.length$gene.length <- gene.length$Gene_End - gene.length$Gene_Start + 1
  gene.length <- aggregate(. ~ Gene_Name, data = gene.length %>% 
                             dplyr::select(matches("Gene_Name|gene.length")), 
                           mean)
  dat.agg <- aggregate(. ~ Gene_Name, data = data.annot %>% 
                           dplyr::select(matches("Gene_Name|GSM")), mean)
  dat.agg <- inner_join(x = dat.agg, y = gene.length, by = "Gene_Name")

  ## Gabel Plot -- Whole
  WT <- c(2,3,4)
  KO <- c(5,6,7)
  dat.agg <- logofMeans.between.A.B(dat = dat.agg, WT, KO)
  p1.all <- gabels.plot(dat.agg[,c("logFC.crude","gene.length")], 
                        comp.between = "(KO/WT)")
  p1.all$plot + coord_cartesian(ylim = c(-0.15,0.15))
  dat.agg$log2.KO.WT1 <- apply(dat.agg[,c(2:5)], 1, 
                               function(r){log2((r[2]+1)/(r[1]+1))})
  dat.agg$log2.KO.WT2 <- apply(dat.agg[,c(3:6)], 1, 
                               function(r){log2((r[2]+1)/(r[1]+1))})
  dat.agg$log2.KO.WT3 <- apply(dat.agg[,c(4:7)], 1, 
                               function(r){log2((r[2]+1)/(r[1]+1))})
  g1 <- gabels.plot(dat.agg[,c("log2.KO.WT1","gene.length")], 
                    comp.between = "(KO/WT)")
  g2 <- gabels.plot(dat.agg[,c("log2.KO.WT2","gene.length")], 
                    comp.between = "(KO/WT)")
  g3 <- gabels.plot(dat.agg[,c("log2.KO.WT3","gene.length")], 
                    comp.between = "(KO/WT)")
  # RTT.Deng <- list(dat.annot = dat.agg, mds.plot = mds.plot, 
  #                  mecp2.plot = mecp2.plot, eset = eset,
  #                  gabel.plots = list(yrs.2 = g1, yrs.5 = g2, yrs.10 = g3, 
  #                                     all = p1.all))
  # save(file = "../results/RTT_Deng.RData",RTT.Deng)
  return(RTT.Deng)
}

## Callosal Neurons RTT Dataset
kishi.callosalneurons.dataanalysis <- function(){
  cat("Running Kishi's Callosal Projection Neurons Mecp2 Dataset\n") 
  experiment <- new("MIAME", name="Kishi Callosal Projection Neurons 
                    RTT Dataset", lab="Jeffrey D. Macklis Lab, Harvard", 
                    title="Reduction of aberrant NF-κB signalling ameliorates 
                    Rett syndrome phenotypes in Mecp2-null mice")
  annotation <- "GPL1261/Affymetrix Mouse Genome 430 2.0 Array"
  load(file = "../dat-infofiles/mm9.GPL1261.RData")
  library("pd.mouse430.2")
  expt <- "GSE50225_RAW"
  group <- c("G1","G1","G1","G2","G2","G2")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  genotypes <- factor(c(rep("WT", 3), rep("KO",3)),levels = c("KO", "WT"))
  mice.strain <- c(rep("C57BL",6))
  sample.info <- factor(c(rep("WT_C57BL", 3), rep("KO_C57BL", 3)),
                        levels = c("KO_C57BL", "WT_C57BL"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(
                                    sample.name = colnames(dat$norm.data),
                                    sample.type = sample.info,
                                    mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData,
                        experimentData = experiment, annotation = annotation)
  mecp2.id = c("1460246_at","MeCP2")
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, samples = sample.info)
  dat1 <- 2^dat$norm.data
  mds.plot <- MDSplot(data = dat1, genotypes, mice.strain)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
  #RTT.callosal.neurons <- list(dat.annot = dat.annot, mds.plot = mds.plot, 
  #                             mecp2.plot = mecp2.plot, eset = eset)
  #save(file = "../results/RTT_callosal.neurons.RData",RTT.callosal.neurons)
  return(list(dat.annot = dat.annot, mds.plot = mds.plot, 
              mecp2.plot = mecp2.plot, eset = eset))
}

## Ben-Shachar Dataset
ben.cere.dataanalysis <- function(){
  cat("Running Ben-Shachar's Cerebellum Dataset\n") 
  load("../dat-infofiles/mm9.GPL6193.RData")
  experiment <- new("MIAME", name="RTT Cerebellum Dataset", 
                    lab="Huda Zoghi Lab, BCM", 
                    title="Mouse models of MeCP2 disorders share gene 
                    expression changes in the cerebellum and hypothalamus.")
  annotation <- "GPL6193/Affymetrix Mouse Exon 1.0 ST Array"
  expt <- "GSE15574_RAW"
  group <- c("G0","G1","G1","G1","G0","G0","G0","G1","G1","G0","G2","G3","G2",
             "G3","G2","G3","G2","G3","G3","G2")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  genotypes <- factor(c(rep("Tg3", 5), rep("WT",5), rep("KO", 5), rep("WT",5)),
                      levels = c("KO", "WT", "Tg3"))
  mice.strain <- c(rep("FVB", 10), rep("C57BL",10))
  sample.info <- factor(c(rep("Tg3_FVB", 5), rep("WT_FVB", 5), 
                          rep("KO_C57BL", 5), rep("WT_C57BL", 5)),
                        levels = c("KO_C57BL", "WT_C57BL", "WT_FVB", "Tg3_FVB"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(
                                        sample.name = colnames(dat$norm.data),
                                        sample.type = sample.info,
                                        mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData,
                        experimentData = experiment, 
                        annotation = annotation)
  title <- "GSE15574/GPL6193 samples"
  boxPlot(data = dat$norm.data, title = title, samples = sample.info)
  mecp2.id = c("7017610","MeCP2") 
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, samples = sample.info)
  dat1 <- 2^dat$norm.data
  mds.plot <- MDSplot(data = dat$norm.data, genotypes, mice.strain)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
  #RTT.cerebellum <- list(dat.annot = dat.annot, mds.plot = mds.plot, 
  #                       mecp2.plot = mecp2.plot, eset = eset)
  #save(file = "../results/RTT_cerebellum.RData",RTT.cerebellum)
  return(list(dat.annot = dat.annot, mds.plot = mds.plot, 
              mecp2.plot = mecp2.plot, eset = eset))
}

## Sugino Neuron Dataset
sugino.neuron.dataanalysis <- function(){
  cat("Running Ben-Shachar's Cerebellum Dataset\n") 
  load("../dat-infofiles/mm9.GPL6193.RData")
  experiment <- new("MIAME", name="RTT Cerebellum Dataset", 
                    lab="Huda Zoghi Lab, BCM", 
                    title="Mouse models of MeCP2 disorders share gene expression
                    changes in the cerebellum and hypothalamus.")
  annotation <- "GPL6193/Affymetrix Mouse Exon 1.0 ST Array"
  expt <- "GSE15574_RAW"
  group <- c("G0","G1","G1","G1","G0","G0","G0","G1","G1","G0","G2","G3","G2",
             "G3","G2","G3","G2","G3","G3","G2")
  dat <- read.geo.function(expt = expt, group = group)
  group <- dat$sort.group
  genotypes <- factor(c(rep("Tg3", 5), rep("WT",5), rep("KO", 5), rep("WT",5)),
                      levels = c("KO", "WT", "Tg3"))
  mice.strain <- c(rep("FVB", 10), rep("C57BL",10))
  sample.info <- factor(c(rep("Tg3_FVB", 5), rep("WT_FVB", 5), 
                          rep("KO_C57BL", 5), rep("WT_C57BL", 5)),
                        levels = c("KO_C57BL", "WT_C57BL", "WT_FVB", "Tg3_FVB"))
  phenoData <- new("AnnotatedDataFrame", data = data.frame(
                                        sample.name = colnames(dat$norm.data),
                                        sample.type = sample.info,
                                        mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat$norm.data)
  eset <- ExpressionSet(assayData = dat$norm.data, phenoData = phenoData,
                        experimentData = experiment, annotation = annotation)
  title <- "GSE15574/GPL6193 samples"
  boxPlot(data = dat$norm.data, title = title, samples = sample.info)
  mecp2.id = c("7017610","MeCP2") 
  mecp2.plot <- GeneExpLevels(data = dat$norm.data, gene.id = mecp2.id, 
                              genotypes = genotypes, samples = sample.info)
  dat1 <- 2^dat$norm.data
  mds.plot <- MDSplot(data = dat$norm.data, genotypes, mice.strain)
  dat.annot <- merge.dat.annot(exprs.dat = dat1, annot.mat = featureData)
  #RTT.cerebellum <- list(dat.annot = dat.annot, mds.plot = mds.plot, 
  #                       mecp2.plot = mecp2.plot, eset = eset)
  #save(file = "../results/RTT_cerebellum.RData",RTT.cerebellum)
  return(list(dat.annot = dat.annot, mds.plot = mds.plot, 
              mecp2.plot = mecp2.plot, eset = eset))
}