#############################################################################
# @Author: Ayush T Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 19th April 2018
#
# Program is used for:
# 1. Analysis of RTT Human Negative Dataset -- Figure 2
#############################################################################

figure2 <- function(){
    load("../results/RTT_Bill.RData")
    cat("\n\n Fig. 2(A) -- Lowry's Human RTT in vitro Dataset \n\n")
    cat("\n\n iPSC Dataset (RTT/WT) \n\n")
    print(RTT.bill$iPSCplot)
    ggsave(filename = "../../plots/Figure2A_1.png",plot = RTT.bill$iPSCplot,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    cat("\n\n NPC Dataset (RTT/WT) \n\n")
    print(RTT.bill$NPCplot)
    ggsave(filename = "../../plots/Figure2A_2.png",plot = RTT.bill$NPCplot,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    cat("\n\n Neuron Dataset (RTT/WT) \n\n")
    print(RTT.bill$Neuron.plot)
    ggsave(filename = "../../plots/Figure2A_3.png",plot = RTT.bill$Neuron.plot,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    cat("\n\n Fig. 2(B) -- Deng's Dataset (RTT/WT) \n\n")
    deng.human.analysis()
    load("../results/RTT_FT.human.RData")
    cat("\n\n Fig. 2(C) -- Lin's Dataset (RTT/WT) \n\n")
    cat("\n\n Frontal Cortex \n\n")
    print(RTT.FT.human$FC + coord_cartesian(ylim = c(-0.15,0.15)))
    cat("\n\n Frontal Cortex \n\n")
    print(RTT.FT.human$TC + coord_cartesian(ylim = c(-0.15,0.15)))
    cat("\n\n Fig. 2(D) -- Lowry's Human RTT Dataset (RTT/WT) \n\n")
    cat("\n\n Frontal Cortex (18 years; Female) \n\n")
    print(RTT.bill$FC1)
    ggsave(filename = "../../plots/Figure2D_1.png",plot = RTT.bill$FC1,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
    cat("\n\n Frontal Cortex (1 year; Male) \n\n")
    print(RTT.bill$FC2)
    ggsave(filename = "../../plots/Figure2D_2.png",plot = RTT.bill$FC2,
           width = 6.50, height = 6.50, dpi = 600, type = "cairo")
}