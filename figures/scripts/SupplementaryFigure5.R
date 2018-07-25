library('splatter')
library('scran')
library('scDD')
library('pscl')
library('ggpubr')
library(ggplot2)

load("../data/SupplementaryFigure5_lunsim.Rdata")
load("../data/SupplementaryFigure5_lun2sim.Rdata")
load("../data/SupplementaryFigure5_simplesim.Rdata")
load("../data/SupplementaryFigure5_Kallisto.Rdata")
load("../data/SupplementaryFigure5_RSEM.Rdata")

comparison<-compareSCEs(list(Simple=simplesim, Lun=lunsim, Lun2=lun2sim, RSEM=RSEM, real=Kallisto))

ggarrange(comparison$Plots$Means, comparison$Plots$Variances, comparison$Plots$MeanVar, comparison$Plots$LibrarySizes, comparison$Plots$ZerosGene, comparison$Plots$ZerosCell, comparison$Plots$MeanZeros, ncol=2, nrow=4)
ggsave("../pngs/SupplementaryFigure5.png", last_plot(), width = 8.64 * 1.5, height = 6.03 * 1.5)
