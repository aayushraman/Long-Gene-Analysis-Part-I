## Illumina annotation file
rm(list =ls())
R1 <- read.table(
    file = "../dat-infofiles/annot.files/IlluminaHumanHT-12.V4.R1.txt",
    sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings = FALSE, 
    fill = TRUE)
R2 <- read.table(
    file = "../dat-infofiles/annot.files/IlluminaHumanHT-12.V4.R2.txt",
    sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings = FALSE, 
    fill = TRUE)
colnames(R1)
colnames(R2)

## updated probes
sum(R1$V14 %in% R2$V14)
idx <- which(R1$V14 %in% R2$V14)
R1 <- R1[-idx, ]
R3 <- rbind(R1,R2)
colnames(R3) <- R3[1,]
R3 <- R3[-1,]

## merge with RefSeq annot
R3$RefSeq_ID <- gsub(pattern = "\\.(.*)", replacement = "",x = R3$RefSeq_ID)
hg19.annot <- read.table(
    file = "../dat-infofiles/annot.files/RefSeq.hg19.Gene.TransciptLength.txt", 
    header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
featureData <- inner_join(hg19.annot, R3, by = c("RefSeq.ID" = "RefSeq_ID"))
colnames(featureData)
featureData$Chromosome <- gsub(pattern = "(.*)", replacement = "chr\\1",
                               x = featureData$Chromosome)
idx <- which(featureData$Chr_Name == featureData$Chromosome)
featureData <- featureData[idx,]
featureData <- featureData[,c(1:9,13,14,22)]
colnames(featureData)
write.table(x = featureData, 
        file = "../dat-infofiles/annot.files/IlluminaHumanHT-12.V4.hg19.txt", 
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)