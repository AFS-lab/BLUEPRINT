library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)

#######################################################################
# FUNCTIONS

#function removes trailing strings from tool names
remove_trail<-function(ggplot_results){
  ggplot_results$Tool<-sub("_cor", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_F1", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_nmrse", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_precision", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_recall", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("Salmon_align", "Salmon\nAlign", ggplot_results$Tool)
  ggplot_results$Tool<-sub("Salmon_quasi", "Salmon\nQuasi", ggplot_results$Tool)
  ggplot_results$Tool<-sub("Salmon_SMEM", "Salmon\nSMEM", ggplot_results$Tool)
  return(ggplot_results)
}

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#function to create ggplot object
make_ggplot<-function(df, title, ylabel){
  p<-ggplot(data=df, aes(x=Tool, y=Value, colour=Experiment)) + geom_point( position=position_jitter(width = .2), stat = "identity") + facet_grid(~Tool, scales= "free_x",space = "free_x")
  p<-p + theme(legend.position = 'none', axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(), text = element_text(size=14), strip.text.x = element_text(size=14))
  p<-p + ggtitle(title) + ylab(ylabel) + scale_colour_manual(values=cbbPalette)
  return(p)
}

#################################

#TO DO

#READ IN BULK DATA
B_bulk<-read.table("../data/Figure4.txt")

#delete Sample_name column and create an Experiment column
B_bulk<-cbind(B_bulk[,1:2], Value=B_bulk[,4], Experiment="bulk")

#Read in Blueprint
B_single<-read.table("../data/Figure1.txt")

#Only keep RSEM Blueprint simulations
B_single<-B_single[B_single$Simulation=="RSEMsim",]

#Delete Simulation column and create an Experiment column
B_single<-cbind(B_single[,2:4], Experiment="single")

#Fuse bulk and single cell dfs
B_df<-rbind(B_bulk, B_single)

#Remove trailing strings after tool names
B_df<-remove_trail(B_df)

#Figuren style Blueprint plots
B_spear<-B_df[B_df$Statistic=="spearmans",]
B_nrmse<-B_df[B_df$Statistic=="nrmse",]
B_precision<-B_df[B_df$Statistic=="precision",]
B_recall<-B_df[B_df$Statistic=="recall",]
B_F1<-B_df[B_df$Statistic=="F1",]


spearmans<-make_ggplot(B_spear, "Spearman's Rho", "Spearman's Rho")
nrmse<-make_ggplot(B_nrmse, "NRMSE", "NRMSE")
precision<-make_ggplot(B_precision, "Precision", "Precision")
recall<-make_ggplot(B_recall, "Recall", "Recall")
F1<-make_ggplot(B_F1, "F1", "F1")

ggarrange(F1, ggarrange(precision, recall, nrow=2), spearmans,nrmse,                                         
          nrow = 2,
          ncol = 2,
          labels = c("A","","B", "C")                                 
) 

ggsave("../pdfs/Figure4.pdf", plot = last_plot(), width=170 *2, units=c("mm") )
ggsave("../pngs/Figure4.png", plot = last_plot(), width=170 *2, units=c("mm") )