## Extracting genomic coordinates
genomicCoord <- function(x){
  if(length(grep("NCBIM37", unlist(x)) > 0)){ ## NCBIM37 for mm9, GRCh38 for hg19, GRCh37 
    ind = grep("NCBIM37", unlist(x))[[1]]  ## NCBIM37 for mm9, GRCh38 for hg19, GRCh37
    gc = x[ind]
    gc = gsub(pattern = "(ensembl_havana_transcript|cdna|havana|ensembl|ensembl_lincrna|ncrna|havana_ig_gene):.*:NCBIM37:([0-9]+|MT|X|Y):([0-9]+):([0-9]+):(1|-1) .*", replacement = "\\2:\\3-\\4", x = gc) ## "(cdna|havana|ensembl_lincrna|ncrna):.*:NCBIM37:([0-9]+|MT):([0-9]+):([0-9]+):(1|-1) .*" for mm9
  }
  else{
    gc = 0
  }
  return(gc)
} 

## Extract probe ids, entrez symbols, and entrez ids
mappingProbe2Gene <- function(featureData){
  # featureData <- read.table("../dat-infofiles/HuGene-2_0-st-v1.na34.hg19.transcript.csv/HuGene-2_0-st-v1.na34.hg19.transcript--withoutheader.txt",
  #                           stringsAsFactors = FALSE, na.strings = FALSE, sep ="\t", quote = "", header = TRUE, fill = TRUE)

  dim(featureData)
  rownames(featureData) <- featureData$probeset_id ## featureData$Probe.Set.ID
  temp1 <- strsplit(featureData[,"mrna_assignment"], "//") ## check gene_assignment mrna_assignment
  temp2 <- unlist(lapply(temp1, function(x) genomicCoord(x)))
  featureData$genomic.coord = temp2
  ind <- featureData$genomic.coord != 0
  featureData$chr <- ifelse(ind == TRUE, sub(pattern = " ([0-9]+|X|Y|MT):([0-9]+)-([0-9]+)", replacement = "chr\\1", x = featureData$genomic.coord, perl = TRUE), NA)
  featureData$gene.start <- ifelse(ind == TRUE, sub(pattern = " ([0-9]+|X|Y|MT):([0-9]+)-([0-9]+)", replacement = "\\2", x = featureData$genomic.coord, perl = TRUE), NA)
  featureData$gene.end <- ifelse(ind == TRUE, sub(pattern = " ([0-9]+|X|Y|MT):([0-9]+)-([0-9]+)", replacement = "\\3", x = featureData$genomic.coord, perl = TRUE), NA)
  featureData$gene.length <- as.numeric(featureData$gene.end) - as.numeric(featureData$gene.start) + 1
  featureData <- featureData[!is.na(featureData$gene.length), ]
  
  ## Extracting the geneName
  temp1 <- strsplit(featureData[,"gene_assignment"], "//")
  featureData$gene.name <- unlist(lapply(temp1, function(x) x[2]))
  featureData$gene.name <- gsub(pattern = " (.*) ", replacement = "\\1", x = featureData$gene.name)
  featureData <- featureData[!is.na(featureData$gene.name), ]
  
  ## removing the probes that may have wrong genomic annotation
  featureData <- featureData[which(featureData$seqname == featureData$chr),]
  #save(featureData, file = "../dat-infofiles/mm9.GPL6193.RData")
  #save(featureData, file = "../dat-infofiles/mm9.GPL6246.RData")
}