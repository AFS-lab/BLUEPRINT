library(ggplot2)
library(ggpubr)
library(reshape2)
library(hydroGOF)
library(scater)

####################################
# FUNCTIONS

#Function that returns legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

#Function for making graphs of results
plot_data<-function( df, title, ylabel, xlabel, ylimits) {
  spearmans_graph<-ggplot(data=df, aes(x=Tool, y=Value, group=Simulation, fill=factor(Simulation)))  + geom_bar(position="dodge", stat = "summary", fun.y = "mean")
  spearmans_graph<- spearmans_graph + scale_x_discrete(labels=c("RSEM", "Salmon Alignment", "Salmon Quasi", "Salmon SMEM", "Sailfish", "eXpress", "Kallisto"))
  spearmans_graph<-spearmans_graph + coord_cartesian(ylim=ylimits) + labs(fill="Simulation Method") + scale_fill_manual(labels = c("Polyester with 3' bias", "RSEM", "Polyester without 3' bias"), values=c("#E69F00", "#56B4E9", "#009E73"))
  spearmans_graph<-spearmans_graph + theme(text=element_text(size=12), legend.position = 'none',axis.text.x=element_text( angle=30,vjust=.8, hjust=0.8)) + ylab(ylabel) + ggtitle(title) + xlab(xlabel)
  return(spearmans_graph)
}


##############################################
# MAKE FIGURE 1

#Read in performance statistics
ggplot_results<-read.table("../data/Figure1.txt")

#Set levels for plotting
ggplot_results$Simulation<-factor(ggplot_results$Simulation, levels=c("RSEMsim", "bias","unbias"))
ggplot_results$Tool<-sub("_cor", "", ggplot_results$Tool)
ggplot_results$Tool<-sub("_F1", "", ggplot_results$Tool)
ggplot_results$Tool<-sub("_nmrse", "", ggplot_results$Tool)
ggplot_results$Tool<-sub("_precision", "", ggplot_results$Tool)
ggplot_results$Tool<-sub("_recall", "", ggplot_results$Tool)
ggplot_results$Tool<-factor(ggplot_results$Tool, levels=c("RSEM", "Salmon_align","Salmon_quasi", "Salmon_SMEM", "Sailfish", "eXpress", "Kallisto"))

#Split ggplot_results by statistic
spearmans<-ggplot_results[ggplot_results$Statistic=="spearmans",]
nrmse<-ggplot_results[ggplot_results$Statistic=="nrmse",]
precision<-ggplot_results[ggplot_results$Statistic=="precision",]
recall<-ggplot_results[ggplot_results$Statistic=="recall",]
F1<-ggplot_results[ggplot_results$Statistic=="F1",]

#make legend
spearmans_graph<-ggplot(data=spearmans, aes(x=Tool, y=Value, group=Simulation, fill=factor(Simulation)))  + geom_bar(position="dodge", stat = "summary", fun.y = "mean")
spearmans_graph<- spearmans_graph + scale_x_discrete(labels=c("RSEM", "Salmon Alignment", "Salmon Quasi", "Salmon SMEM", "Sailfish", "eXpress", "Kallisto"))
spearmans_graph<-spearmans_graph + labs(fill="Simulation Method") + scale_fill_manual(labels =c( "RSEM", "Polyester (3' bias)", "Polyester (no 3' bias)"), values=c("#E69F00", "#56B4E9", "#009E73")) + theme(text=element_text(size=12))
leg<-g_legend(spearmans_graph)

#make graphs
spearmans_graph<-plot_data(spearmans, "Spearman's Rho", "Spearman's Rho", "", c(0,1))
nrmse_graph<-plot_data(nrmse, "NRMSE", "NRMSE", "", c(0,70))
precision_graph<-plot_data(precision,"Precision", "Precision", "", c(0,1))
recall_graph<-plot_data(recall, "Recall", "Recall", "", c(0,1))
F1_graph<-plot_data(F1, "F1", "F1", "", c(0,1))

#arrange in figure
ggarrange(F1_graph,ggarrange(precision_graph,recall_graph, nrow=2),leg, spearmans_graph,nrmse_graph,                                        
          nrow = 2,
          ncol = 3,
          labels = c("A","","","B","C",""),
          widths=c(2,2,2)
) 

ggsave("../pdfs/Figure1.pdf", plot=last_plot(), height= 225, width=170, units=c("mm"))
ggsave("../pngs/Figure1.png", last_plot(), height= 225, width=170, units=c("mm"))
