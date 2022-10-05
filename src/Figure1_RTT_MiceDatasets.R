#############################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 19th April 2018
#
# Program is used for:
# 1. Analysis of RTT Microarray Dataset -- Figure 1
#############################################################################

figure1 <- function(){

    # Topotecan -- King            
    figure1_topotecan()
    
    # RTT Amygdala dataset
    load("../results/RTT_amygdala.RData")
    
    ## Run the function to reproduce the results
    # RTT.amygdala <- samaco.amygdala.dataanalysis()
  
    ## overlay gabel plot for WT/WT and KO/WT C57BL 
    log2FC.length.WT <- log2FCwithingenotypes(RTT.amygdala$dat.annot[,c(1:5,20)])
    log2FC.length.KO <- logofMeans.between.A.B(dat = RTT.amygdala$dat.annot, 
                                               A.samples = c(2:5), 
                                               B.samples = c(6:9))
    log2FC.length <- inner_join(x=log2FC.length.WT[,c("gene.name","comp.mat2")],
                                y = log2FC.length.KO[,c("gene.name",
                                                        "logFC.crude",
                                                        "gene.length")], 
                                by = "gene.name")
    p1 <- overlay.gabels.plot(log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)",
                              comp.between2 = "(KO/WT)") 
    cat("\n\n Fig. 1(B) -- Amygdala (KO/WT)\n\n")
    print(plot_grid(p1$plot1 + coord_cartesian(ylim = c(-0.25,0.15)), 
                      p1$plot2 + coord_cartesian(ylim = c(0,10)),
                      ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
    p1 <- plot_grid(p1$plot1 + coord_cartesian(ylim = c(-0.25,0.15)), 
                    p1$plot2 + coord_cartesian(ylim = c(0,10)),
                    ncol = 1, align = 'v') + 
        theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
    ggsave(filename = "../../plots/Figure1B_1.png",plot = p1,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length, log2FC.length.WT, log2FC.length.KO)
  
    ## overlay gabel plot WT/WT and Tg/WT for FVB
    log2FC.length.WT <- log2FCwithingenotypes(RTT.amygdala$dat.annot[,c(1,15:18,20)])
    log2FC.length.Tg <- logofMeans.between.A.B(dat = RTT.amygdala$dat.annot,
                                               A.samples = c(15:19), 
                                               B.samples = c(10:14))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name",
                                                        "comp.mat1")], 
                                y = log2FC.length.Tg[,c("gene.name",
                                                        "logFC.crude",
                                                        "gene.length")], 
                                by = "gene.name")
    p2 <- overlay.gabels.plot(log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)",
                              comp.between2 = "(Tg/WT)") 
    cat("\n\n Fig. 1(C) -- Amygdala (Tg/WT)\n\n")
    print(plot_grid(p2$plot1 + coord_cartesian(ylim = c(-0.15,0.25)), p2$plot2 + 
                          coord_cartesian(ylim = c(0,10)), ncol = 1,
                    align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), 
                                                            "cm")))
    p2 <- plot_grid(p2$plot1 + coord_cartesian(ylim = c(-0.15,0.25)), p2$plot2 + 
                        coord_cartesian(ylim = c(0,10)), ncol = 1,
                    align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), 
                                                            "cm"))
    ggsave(filename = "../../plots/Figure1C_1.png",plot = p2,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length.WT, log2FC.length.Tg, log2FC.length)
  
    # RTT Cerebellum dataset
    load("../results/RTT_cerebellum.RData")
    
    ## Run the function to reproduce the results
    #RTT.cerebellum <- ben.cere.dataanalysis()
  
    ## overlay gabel plot for WT/WT and KO/WT for C57BL
    log2FC.length.WT <- log2FCwithingenotypes(RTT.cerebellum$dat.annot[,c(1,18:22)])
    log2FC.length.KO <- logofMeans.between.A.B(dat = RTT.cerebellum$dat.annot, 
                                    A.samples = c(17:21), B.samples = c(12:16))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat1")], 
                              y = log2FC.length.KO[,c("gene.name","logFC.crude",
                                            "gene.length")], by = "gene.name")
    p3 <- overlay.gabels.plot(log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)",
                              comp.between2 = "(KO/WT)")
    cat("\n\n Fig. 1(B) -- Cerebellum (KO/WT)\n\n")
    print(plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p3$plot2 + 
                        coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
    p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p3$plot2 + 
                        coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v')+ 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
    ggsave(filename = "../../plots/Figure1B_2.png",plot = p3,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length, log2FC.length.WT, log2FC.length.KO)
  
    ## overlay gabel plot for WT/WT and KO/WT FVB samples
    log2FC.length.WT <- log2FCwithingenotypes(RTT.cerebellum$dat.annot[,c(1,7,9:11,22)])
    log2FC.length.Tg <- logofMeans.between.A.B(dat = RTT.cerebellum$dat.annot, 
                                        A.samples = c(7:11), B.samples = c(2:6))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat3")], 
                              y = log2FC.length.Tg[,c("gene.name","logFC.crude", 
                                                      "gene.length")], 
                              by = "gene.name")
    p4 <- overlay.gabels.plot(log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", 
                              comp.between2 = "(Tg/WT)")
    cat("\n\n Fig. 1(C) -- Cerebellum (Tg/WT)\n\n")
    print(plot_grid(p4$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p4$plot2 + 
                        coord_cartesian(ylim = c(0,10)),ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
    p4 <- plot_grid(p4$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p4$plot2 + 
                        coord_cartesian(ylim = c(0,10)),ncol = 1, align = 'v') + 
        theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
    ggsave(filename = "../../plots/Figure1C_2.png",plot = p4,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length, log2FC.length.WT, log2FC.length.Tg)
  
    # RTT Hypothalamus dataset
    load("../results/RTT_hypo.RData")
    ## Run the function to reproduce the results
    #  RTT.hypo <- maria.hypothalamus.dataanalysis()
  
    ## overlay gabel plot for WT/WT and KO/WT C57BL
    log2FC.length.WT <- log2FCwithingenotypes(dat = RTT.hypo$dat.annot[,c(1,6:9,18)])
    log2FC.length.KO <- logofMeans.between.A.B(dat = RTT.hypo$dat.annot, 
                                    A.samples = c(6:9), B.samples = c(10:13))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat3")],
                                y = log2FC.length.KO[,c("gene.name",
                                                        "logFC.crude",
                                                        "gene.length")], 
                                by = "gene.name")
    p5 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], comp.between1 = 
                                  "(WT/WT)", comp.between2 = "(KO/WT)")
    cat("\n\n Fig. 1(B) -- Hypothalamus (KO/WT) \n\n")
    print(plot_grid(p5$plot1 + coord_cartesian(ylim = c(-0.15,0.15)),p5$plot2 + 
                    coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
    p5 <- plot_grid(p5$plot1 + coord_cartesian(ylim = c(-0.15,0.15)),p5$plot2 + 
                        coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v')+ 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
    ggsave(filename = "../../plots/Figure1B_3.png",plot = p5,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length, log2FC.length.WT, log2FC.length.KO)
  
    ## overlay gabel plot for WT/WT (GSM281260, GSM281261, GSM281265, 
    ## GSM281266 -- check box plot) and Tg/WT FVB
    
    log2FC.length.WT <- log2FCwithingenotypes(dat = RTT.hypo$dat.annot[,c(1:5,18)])
    log2FC.length.Tg <- logofMeans.between.A.B(dat = RTT.hypo$dat.annot,
                                               A.samples = c(2:5),
                                               B.samples = c(14:17))
    log2FC.length.WT$comp.mat2 <- -log2FC.length.WT$comp.mat2
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat2")], 
                              y = log2FC.length.Tg[,c("gene.name","logFC.crude",
                                                      "gene.length")], 
                              by = "gene.name")
    p6 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], 
                              comp.between1 = "(WT/WT)", 
                              comp.between2 = "(Tg/WT)")
    cat("\n\n Fig. 1(C) -- Hypothalamus (Tg/WT) \n\n")
    print(plot_grid(p6$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p6$plot2 +
                    coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
    p6 <- plot_grid(p6$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p6$plot2 +
                  coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v') + 
        theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
    ggsave(filename = "../../plots/Figure1C_3.png",plot = p6,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length, log2FC.length.WT, log2FC.length.Tg)

    # Johnson et al. Nature Medicine 2017
    load("../results/RTT_Johnson_Neurons.RData")
    cat("\n\n Fig. 1(D) -- R106W Excitatory Neurons Female \n\n")
    print(Johnson_Mecp2[[1]])
    ggsave(filename = "../../plots/Figure1D_1.png",plot = Johnson_Mecp2[[1]],
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    cat("\n\n Fig. 1(D) -- R106W MUT Excitatory Neurons Female \n\n")
    print(Johnson_Mecp2[[2]])
    ggsave(filename = "../../plots/Figure1D_2.png",plot = Johnson_Mecp2[[2]],
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    cat("\n\n Fig. 1(D) -- T158M MUT Excitatory Neurons Female \n\n")
    print(Johnson_Mecp2[[3]])
    ggsave(filename = "../../plots/Figure1D_3.png",plot = Johnson_Mecp2[[3]],
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
}

figureS2 <- function(){
    
    # RTT Straitum dataset
    load("../results/RTT_striatum.RData")
    ## Run the function to reproduce the results
    #RTT.striatum <- zhao.striatum.dataanalysis()
    
    ## overlay gabel plot WT/WT for C57BL
    log2FC.length.WT <- log2FCwithingenotypes(RTT.striatum$dat.annot[,c(1,3:6,22)])
    log2FC.length.KO <- logofMeans.between.A.B(dat = RTT.striatum$dat.annot, 
                                               A.samples = c(2:6), 
                                               B.samples = c(7:11))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat1")], 
                                y = log2FC.length.KO[,c("gene.name",
                                                        "logFC.crude",
                                                        "gene.length")], 
                                by = "gene.name")
    p7 <- overlay.gabels.plot(log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", 
                              comp.between2 = "(KO/WT)")
    cat("\n\n Supplemnetary Fig. 2(A)-- Striatum\n\n")
    print(plot_grid(p7$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p7$plot2 + 
                        coord_cartesian(ylim = c(0,10)), ncol = 1, align = 'v')+ 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
    rm(log2FC.length.WT, log2FC.length.KO, log2FC.length)
    
    # RTT Hippocampus dataset
    load("../results/RTT_hippo.RData")
    ## Run the function to reproduce the results
    #RTT.hippo <- baker.hippocampus.dataanalysis()
    
    ## overlap gabel plot for WT/WT FVBx129 and KO/WT -- 4 weeks
    log2FC.length.WT <- log2FCwithingenotypes(RTT.hippo$dat.annot[,c(1,2:5,34)])
    log2FC.length.KO <- logofMeans.between.A.B(dat = RTT.hippo$dat.annot, 
                                               A.samples = c(2:5), 
                                               B.samples = c(6:9))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat2")], 
                                y = log2FC.length.KO[,c("gene.name", 
                                                        "logFC.crude",
                                                        "gene.length")],
                                by = "gene.name")
    p8 <- overlay.gabels.plot(log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", 
                              comp.between2 = "(KO/WT)")
    cat("\n\n Supplemnetary Fig. 2(B) -- Hippocampus (4 weeks) \n\n")
    print(plot_grid(p8$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p8$plot2 + 
                        coord_cartesian(ylim = c(0,10)), ncol = 1, 
                    align = 'v') + theme(plot.margin = unit(c(1,0.5,0.5,0.5), 
                                                            "cm")))
    rm(log2FC.length, log2FC.length.WT, log2FC.length.KO)
    
    ## overlap plot for WT/WT FVBx129 and KO/WT -- 9 weeks
    log2FC.length.WT <- log2FCwithingenotypes(RTT.hippo$dat.annot[,c(1,18:21,34)])
    log2FC.length.KO <- logofMeans.between.A.B(dat = RTT.hippo$dat.annot, 
                                    A.samples = c(18:21), B.samples = c(22:25))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name",
                                                        "comp.mat3")], 
                                y = log2FC.length.KO[,c("gene.name",
                                                        "logFC.crude",
                                                        "gene.length")], 
                                by = "gene.name")
    p9 <- overlay.gabels.plot(log2FC.length[,c(2:4)], comp.between1 = "(WT/WT)", 
                              comp.between2 = "(KO/WT)")
    cat("\n\n Supplemnetary Fig. 2(C) -- Hippocampus (9 weeks) \n\n")
    print(plot_grid(p9$plot1 + coord_cartesian(ylim = c(-0.15,0.15)), p9$plot2 + 
                        coord_cartesian(ylim = c(0,10)), ncol = 1,align = 'v') + 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))
    rm(log2FC.length, log2FC.length.WT, log2FC.length.KO)

    # Gabel Visual Cortex
    load("../results/RTT_Gabel_VC.RData")
    cat("\n\n Supplemnetary Fig. 2(D) -- Visual Cortex (RNA-seq)\n\n")
    print(plot_grid(gabel_Mecp2$p1$plot1 + coord_cartesian(ylim = c(-0.5,0.5)), 
                    gabel_Mecp2$p1$plot2 + coord_cartesian(ylim = c(0,10)), 
                    ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")))

    # Sugino Neuron Dataset
    load("../results/RTT_Sugino_Neurons.RData")
    cat("\n\n Supplemnetary Fig. 2(E) -- Locus Coeruleus (P22) \n\n")
    print(Sugino_Mecp2[[1]])
    cat("\n\n Supplemnetary Fig. 2(F) -- Locus Coeruleus \n\n")
    print(Sugino_Mecp2[[2]])
    cat("\n\n Supplemnetary Fig. 2(G) -- FS, Motor Cortex \n\n")
    print(Sugino_Mecp2[[3]])
    cat("\n\n Supplemnetary Fig. 2(H) -- Purkinje Cells Cerebellum \n\n")
    print(Sugino_Mecp2[[4]])
    cat("\n\n Supplemnetary Fig. 2(I) -- PN, Motor Cortex \n\n")
    print(Sugino_Mecp2[[5]])

    # Callosal Neurons Dataset
    load("../results/RTT_callosal.neurons.RData")
    cat("\n\n Supplemnetary Fig. 2(J) -- Callosal Projection Neurons\n\n")
    print(RTT.callosal.neurons[[1]])

    # Hypo RNA-seq Dataset
    chen.hypothalamus.dataset()
}

figure1_topotecan <- function(){
    load(file = "../results/King_RNAseq.RData")
    log2FC.length.WT <- log2FCwithingenotypes(dat = dat.dds[,c(1,2:5,18)])
    log2FC.length.KO <- logofMeans.between.A.B(dat = dat.dds, A.samples = c(2:6), 
                                               B.samples = c(7:11))
    log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name",
                                                        "comp.mat3")],
                                y = log2FC.length.KO[,c("gene.name",
                                                        "logFC.crude",
                                                        "gene.length")], 
                                by = "gene.name")
    p1 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], 
                              comp.between1 = "(V/V)", comp.between2 = "(D/V)")
    cat("\n\n Fig. 1(A)--Cultured Cortical Neurons; (RNA-Seq; King et al.)\n\n")
    print(plot_grid(p1$plot1 + coord_cartesian(ylim = c(-2,2)), p1$plot2 + 
                    coord_cartesian(ylim = c(0,40)), ncol = 1, align = 'v') + 
                    theme(plot.margin = unit(c(1,1,1,1), "cm")))
    p1 <- plot_grid(p1$plot1 + coord_cartesian(ylim = c(-2,2)), p1$plot2 + 
                    coord_cartesian(ylim = c(0,40)), ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,1,1,1), "cm"))
    ggsave(filename = "../../plots/Figure1A_1.png",plot = p1,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length, log2FC.length.WT, log2FC.length.KO)
    
    load("../results/King_Array.RData")
    dat.annot$comp.mat3 <- apply(X = dat.annot[,c(1,5:8)], 1, 
                            function(r){log2((mean(as.numeric(r[c(2)]))+1)/
                                                 (as.numeric(r[4])+1))})
    log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, 
                                               A.samples = c(5:7), 
                                               B.samples = c(2:4))
    log2FC.length <- inner_join(x = dat.annot[,c("gene.name","comp.mat3")],
                                y = log2FC.length.KO[,c("gene.name",
                                                        "logFC.crude",
                                                        "gene.length")], 
                                by = "gene.name")
    p2 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)], 
                              comp.between1 = "(V/V)", 
                              comp.between2 = "(D/V)")
    cat("\n\n Fig. 1(A) -- Cultured Cortical Neurons; (Array; King et al.)\n\n") 
    print(plot_grid(p2$plot1 + coord_cartesian(ylim = c(-2,2)), p2$plot2 + 
                coord_cartesian(ylim = c(0,40)), ncol = 1, align = 'v') + 
                theme(plot.margin = unit(c(1,1,1,1), "cm")))
    p2 <- plot_grid(p2$plot1 + coord_cartesian(ylim = c(-2,2)), p2$plot2 + 
                        coord_cartesian(ylim = c(0,40)), ncol = 1, align = 'v')+ 
              theme(plot.margin = unit(c(1,1,1,1), "cm"))
    ggsave(filename = "../../plots/Figure1A_2.png",plot = p2,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    rm(log2FC.length, log2FC.length.KO)
    load("../results/Mabb_RNAseq.RData")
    cat("\n\n Fig. 1(A)--Cultured Cortical Neurons; (RNA-seq; Mabb et al.)\n\n") 
    print(plot_grid(Mabb.seq$p1$plot1 + coord_cartesian(ylim = c(-2,2)),
                    Mabb.seq$p1$plot2 + coord_cartesian(ylim = c(0,80)), 
                                                        ncol = 1, align = 'v') + 
            theme(plot.margin = unit(c(1,1,1,1), "cm")))
    p3 <- plot_grid(Mabb.seq$p1$plot1 + coord_cartesian(ylim = c(-2,2)),
                    Mabb.seq$p1$plot2 + coord_cartesian(ylim = c(0,80)), 
                    ncol = 1, align = 'v') + 
              theme(plot.margin = unit(c(1,1,1,1), "cm"))
    ggsave(filename = "../../plots/Figure1A_3.png",plot = p3,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
}

figureS3 <- function(){
    load("../results/RTT_Johnson_Neurons.RData")
    cat("\n\n Supplemnetary Fig. 3 -- GRO-seq \n\n")
    print(Johnson_Mecp2[[5]])
    cat("\n\n Supplemnetary Fig. 3 -- Whole Cell-seq \n\n")
    print(Johnson_Mecp2[[6]])
}

figureS4 <- function(){    
    load("../results/RTT_Johnson_Neurons.RData")
    cat("\n\n Supplemnetary Fig. 4(A) -- R106W Excitatory Nuclear Male \n\n")
    print(Johnson_Mecp2[[7]])
    cat("\n\n Supplemnetary Fig. 4(B) -- T158M Excitatory Nuclear Male \n\n")
    print(Johnson_Mecp2[[8]])
    cat("\n\n Supplemnetary Fig. 4(C) -- R106W Inhibitory Nuclear Male \n\n")
    print(Johnson_Mecp2[[9]])
    cat("\n\n Supplemnetary Fig. 4(D) -- T158M Inhibitory Nuclear Male \n\n")
    print(Johnson_Mecp2[[10]])
    cat("\n\n Supplemnetary Fig. 4(E)--T158M-WT Excitatory Nuclear Female\n\n")
    print(Johnson_Mecp2[[4]])
}