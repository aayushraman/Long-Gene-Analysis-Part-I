## overlay plot
overlay.gabels.plot <- function(mat, length.type = "Gene", comp.between1 = "", 
                                comp.between2 = ""){
  p1 <- overlay.moving.average.function(dat = mat, bin.size = 200,
                                        shift.size = 40, length.type,
                                        comp.between1, comp.between2)
  return(p1)
}  
overlay.moving.average.function <- function(dat, bin.size, shift.size, 
                                            length.type, comp.between1, 
                                            comp.between2){
    dat <- dat[order(dat$gene.length),]
    dat$gene.length <- dat$gene.length/1000
  
    # data frame for storing the values and 
    # calculating the number of bins
    mean.points <- data.frame()
    num.bins <- round((dim(dat)[1]-bin.size)/shift.size)+1
  
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
        mat1 <- dat[start:end, 1]
        mat2 <- dat[start:end, 2]
        mat.mean1 <- mean(mat1)
        mat.mean2 <- mean(mat2)
        mat.sd.1 <- sd(mat1)
        mat.sd.2 <- sd(mat2)
        mat.length <- mean(dat[start:end, 3])
        bin.width <- end-start+1

    ## mat means
    pval <- t.test2(m1 = mat.mean1, m2 = mat.mean2, s1 = mat.sd.1, 
                    s2 = mat.sd.2, n1 = bin.width)
    pval.log10 <- -log10(pval)
    mat.mean <- data.frame(mat.mean1, mat.mean2, mat.length, mat.sd.1, 
                           mat.sd.2, pval, pval.log10)
    mean.points <- rbind(mean.points, mat.mean)
    
    ## end exceeds total number of genes
    if (end == dim(dat)[1]) break;
  }
  ## colors used for moving average plots
  col1 <- comp.between1
  col2 <- comp.between2

  ## overlay line plot
  mean.points$fdr <- p.adjust(p = mean.points$pval, method = "fdr")
  ind <- mean.points$mat.length >=1 & mean.points$mat.length <=1000
  mean.points = mean.points[ind, ]
  plot1 <- ggplot(data = mean.points, aes(x = mat.length)) + 
      geom_line(aes(y = mean.points$mat.mean1, color = col1), size=1) + 
      geom_line(aes(y = mean.points$mat.mean2, color = col2), size=1) + 
      ylab(paste("Mean Log2FC")) + 
      scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
      geom_ribbon(aes(ymin=(mat.mean1-(mat.sd.1*0.50)), 
                      ymax=(mat.mean1+(mat.sd.1*0.50)), 
                      x = mat.length, fill = col1), alpha=.25) +
      geom_ribbon(aes(ymin=(mat.mean2-(mat.sd.2*0.50)), 
                      ymax=(mat.mean2+(mat.sd.2*0.50)), x = mat.length, 
                      fill = col2), alpha=.25) +
      theme(legend.position="none", ## legend.text = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 24, face = "bold"),
            axis.text.y = element_text(size = 24, face = "bold", 
                                       color = "black"),
            axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), axis.line.x = element_blank())
  
    ##p-value plot
    gene.type <- ifelse(test = mean.points$fdr < 0.05, yes = "#FF0000",
                        no = "#696969")
    if(sum(mean.points$fdr < 0.05) > 0){
        y.int <- min(mean.points[which(mean.points$fdr < 0.05 & 
                                           !is.infinite(mean.points$fdr)),
                                 "pval.log10"])
        plot2 <- ggplot(data = mean.points, aes(x = mat.length, 
                                                y = pval.log10)) + 
            geom_line(size = 0.4, colour="gray70") + 
            geom_point(size = 2, color = gene.type) + 
            geom_hline(aes(yintercept = y.int), colour="#FF0000", 
                       linetype="dashed", size = 1) + 
            scale_x_continuous(trans = log10_trans(), 
                               breaks = c(0,1,10,100,1000)) +
            xlab(paste("Mean",length.type,"Length in KB")) + 
            ylab(paste("-Log10(pvalue)")) + theme(legend.position="none",
                        axis.title = element_text(size = 24, face = "bold"),
                        axis.text.x = element_text(size = 24, face = "bold", 
                                                   color = "black"),
                        axis.text.y = element_text(size = 24, face = "bold",
                                                   color = "black"))
        }
    else{
        y.int <- ceiling(max(mean.points[,"pval.log10"]))
        plot2 <- ggplot(data = mean.points, aes(x = mat.length, 
                                                y = pval.log10)) + 
            geom_line(size = 0.4, colour="gray70") + 
            geom_point(size = 2, color = gene.type) + 
            scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,
                                                                 1000)) +
            xlab(paste("Mean",length.type,"Length in KB")) + 
            ylab(paste("-Log10(pvalue)")) + 
            geom_hline(aes(yintercept = y.int), colour="#FF0000", 
                       linetype="dashed", size = 1) + 
            theme(legend.position="none", 
                  axis.title = element_text(size = 24, face = "bold"),
                  axis.text.x = element_text(size = 24, face = "bold", 
                                             color = "black"),
                  axis.text.y = element_text(size = 24, face = "bold", 
                                             color = "black"))
    }
    return(list(plot1 = plot1, plot2 = plot2)) ## ,mat1 = mean.points
}

## p-values from 2 sample t-test; code adapted from http://bit.ly/2eqeYyO
t.test2 <- function(m1,m2,s1,s2,n1,n2=n1,m0=0,equal.variance=FALSE){
    if( equal.variance==FALSE ){
        se <- sqrt( (s1^2/n1) + (s2^2/n2))
        df <- ((s1^2/n1 + s2^2/n2)^2)/((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
        }
    else{
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)) 
        df <- n1+n2-2
    }
    t <- (m1-m2-m0)/se 
    pval <- 2*pt(-abs(t),df)
    return(pval) 
}