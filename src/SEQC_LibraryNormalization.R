## Library Size normalization
All.normalizationMethods <- function(counts.table, A.samples, B.samples, C.samples, hg19.annot){
  A.numCounts <- sum(counts.table[,A.samples])
  B.numCounts <- sum(counts.table[,B.samples])
  C.numCounts <- sum(counts.table[,C.samples])
  A.normfactor <- sum(c(as.numeric(A.numCounts), as.numeric(B.numCounts), as.numeric(C.numCounts)))/A.numCounts
  B.normfactor <- sum(c(as.numeric(A.numCounts), as.numeric(B.numCounts), as.numeric(C.numCounts)))/B.numCounts
  C.normfactor <- sum(c(as.numeric(A.numCounts), as.numeric(B.numCounts), as.numeric(C.numCounts)))/C.numCounts
  
  ##counts after lib size normalization and aggregating the transcripts to genes
  counts.lib.norm <- cbind(counts.table[,A.samples]*A.normfactor, 
                           counts.table[,B.samples]*B.normfactor,
                           counts.table[,C.samples]*C.normfactor)
  counts.lib.norm <- inner_join(x = hg19.annot, y = counts.lib.norm %>%
                                add_rownames(), by = c("RefSeq.ID"="rowname"))
  dat.agg <- merge.dat.annot(exprs.dat = counts.lib.norm)

  ## LN normalization -- gl
  ratio.LN.gl <- logofMeans.between.ABC(dat = dat.agg, A.samples, B.samples, C.samples)
  ratio.LN.gl <- data.frame(ratio.LN.gl[, c("logFC.crude","gene.length", "gene.name")])
  Complot.TC.gl <- gabels.plot(mat = ratio.LN.gl, y.axis = "Mean Log2 (ß ratio)")
  
  ## LN normalization -- tl
  ratio.LN.tl <- logofMeans.between.ABC(dat = counts.lib.norm, A.samples,
                                        B.samples, C.samples)
  ratio.LN.tl <- data.frame(ratio.LN.tl[, c("logFC.crude", "transcript.length", 
                                            "RefSeq.ID")])
  Complot.TC.tl <- gabels.plot(mat = data.frame(ratio.LN.tl[, c("logFC.crude",
                    "transcript.length", "RefSeq.ID")]), 
                    length.type = "Transcript", y.axis = "Mean Log2 (ß ratio)")
  rm(dat.agg)
  
  ##egdeR normalization
  edgeR.obj <- DGEList(counts.table)
  edgeR.obj <- calcNormFactors(edgeR.obj)
  counts.edgeRNorm <- data.frame(cpm(edgeR.obj, log=FALSE))
  counts.edgeRNorm <- inner_join(x = hg19.annot, y = counts.edgeRNorm %>%
                                add_rownames(),by = c("RefSeq.ID"="rowname"))
  dat.agg <- merge.dat.annot(exprs.dat = counts.edgeRNorm)
  
  ## edgeR normalization -- gl
  ratio.edgeR.gl <- logofMeans.between.ABC(dat = dat.agg, A.samples, B.samples, C.samples)
  ratio.edgeR.gl <- data.frame(ratio.edgeR.gl[, c("logFC.crude","gene.length", "gene.name")])
  Complot.TMM.gl <- gabels.plot(mat = ratio.edgeR.gl, y.axis = "Mean Log2 (ß ratio)")
  
  ## edgeR normalization -- tl
  ratio.edgeR.tl <- logofMeans.between.ABC(dat = counts.edgeRNorm, A.samples, 
                                           B.samples, C.samples)
  Complot.TMM.tl <- gabels.plot(mat=data.frame(ratio.edgeR.tl[,c("logFC.crude",
                        "transcript.length", "RefSeq.ID")]), 
                        length.type = "Transcript",y.axis="Mean Log2 (ß ratio)")
  plots <- list(plotA1 = Complot.TC.gl, plotA2 = Complot.TC.tl,
                plotB1 = Complot.TMM.gl, plotB2 = Complot.TMM.tl)
  return(plots)
}