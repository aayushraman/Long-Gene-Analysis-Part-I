Voineagu.frontal.temporal.dataset <- function(){
  
  cat("Running Voineagu I, Lin P's Mecp2 Dataset -- Human Dataset -- GSE75303\n") 
  experiment <- new("MIAME", name="RTT Fronatal/Temporal Cortex Dataset", 
            lab = "I. Voineagu Lab, UNSW", 
            title="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75303")
  annotation <- "GPL10558/Illumina HumanHT-12 V4.0 expression beadchip"

  ## reading the cel files
  dat.norm <- read.table("../dat/GEO/GSE75303_normalized.txt", header = TRUE, 
                         sep = "\t", quote = "", stringsAsFactors = FALSE)
  colnames(dat.norm) <- c("ID_REF","Array_Address_Id","GSM1949097","GSM1949098", 
                         "GSM1949099","GSM1949100", "GSM1949101", "GSM1949102", 
                         "GSM1949103", "GSM1949104", "GSM1949105", "GSM1949106",
                         "GSM1949107", "GSM1949108")
  ind <- which(!is.na(dat.norm$ID_REF))
  dat.norm <- dat.norm[ind,]
  rownames(dat.norm) <- dat.norm$ID_REF
  dat.norm <- as.matrix(dat.norm[,c(3:14)])
  group <- c("G0","G0","G1","G1","G0","G1","G2","G3","G2","G3","G2","G3")
  dat.norm <- dat.norm[ , order(group)]
  group <- group[order(group)]

  ## for the plots
  genotypes <- factor(c(rep("WT",6),rep("RTT",6)),levels = c("RTT", "WT"))
  mice.strain <- c(rep("Frontal Cortex",3), rep("Temporal Cortex",3), 
                   rep("Frontal Cortex",3), rep("Temporal Cortex",3))
  sample.info <- factor(c(rep("FC_WT",3),rep("TC_WT",3),rep("FC_RTT",3),
                          rep("TC_RTT",3)), 
                          levels = c("FC_WT","TC_WT","FC_RTT","TC_RTT"))
  phenoData <- new("AnnotatedDataFrame", data = 
                       data.frame(sample.name = colnames(dat.norm),
                        sample.type = sample.info, mice.strain = mice.strain))
  rownames(phenoData) <- colnames(dat.norm)
  eset <- ExpressionSet(assayData = dat.norm, phenoData = phenoData, 
                        experimentData = experiment, annotation = annotation)
  ## Mecp2 levels
  mecp2.id <- c("ILMN_1702715","MeCP2")
  mecp2.plot1 <- GeneExpLevels(data = dat.norm, gene.id = mecp2.id, 
                               genotypes = genotypes, samples = sample.info)
  mecp2.id <- c("ILMN_3310740","MeCP2")
  mecp2.plot2 <- GeneExpLevels(data = dat.norm, gene.id = mecp2.id, 
                               genotypes = genotypes, samples = sample.info)
  mecp2.id <- c("ILMN_1824898","MeCP2")
  mecp2.plot3 <- GeneExpLevels(data = dat.norm, gene.id = mecp2.id, 
                               genotypes = genotypes, samples = sample.info)
  
  ## PCA plot
  dat1 <- 2^(dat.norm)
  #mds.plot <- MDSplot(data = dat1, genotypes = genotypes, 
  #                    conditions = mice.strain)

  ## merging anot file without gene length
  featureData <- read.table(
      file = "../dat-infofiles/annot.files/IlluminaHumanHT-12.V4.hg19.txt",
                            sep = "\t",header = TRUE,quote = "", fill = TRUE, 
                            stringsAsFactors = FALSE, na.strings = FALSE)
  dat.annot <- merge.dat.annot.illu(exprs.dat = dat1, annot.mat = featureData)
  
  ## KO/WT Frontal Cortex
  log2FC.length <- logofMeans.between.A.B(dat = dat.annot, A.samples = c(2:4), 
                                          B.samples = c(8:10))
  p1 <- gabels.plot(log2FC.length[,c("logFC.crude","gene.length")], 
                        comp.between = "(RTT/WT)")
  log2FC.length <- logofMeans.between.A.B(dat = dat.annot, A.samples = c(5:7), 
                                          B.samples = c(11:13))
  p2 <- gabels.plot(log2FC.length[,c("logFC.crude","gene.length")], 
                    comp.between = "(RTT/WT)") 
  RTT.FT.human <- list(FC = p1, TC = p2)
  save(RTT.FT.human, file = "../results/RTT_FT.human.RData")
  return(RTT.FT.human)
}