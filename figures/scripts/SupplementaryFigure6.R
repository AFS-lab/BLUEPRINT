library("ggplot2")
library("ggpubr")

#######################
#FUNCTIONS

plot_data<-function( df, title, ylabel, xlabel, ylimits) {
  spearmans_graph<-ggplot(data=df, aes(x=Tool, y=Value, group=Simulation, fill=factor(Simulation)))  + geom_bar(position="dodge", stat = "summary", fun.y = "mean")
  spearmans_graph<- spearmans_graph + scale_x_discrete(labels=c("RSEM", "Salmon Alignment Mode", "Salmon Quasi Mode", "Salmon SMEM Mode", "Sailfish", "eXpress", "Kallisto"))
  spearmans_graph<-spearmans_graph + coord_cartesian(ylim=ylimits) + labs(fill="Simulation Method") + scale_fill_manual(labels = c("Polyester with 3' bias", "RSEM", "Polyester without 3' bias"), values=c("#E69F00", "#56B4E9", "#009E73"))
  spearmans_graph<-spearmans_graph + theme(text=element_text(size=16), legend.position = 'none',axis.text.x=element_text( angle=30,vjust=.8, hjust=0.8),plot.margin = unit(c(1,1,1,1), "lines"))+ ylab(ylabel) + ggtitle(title) + xlab(xlabel)
  return(spearmans_graph)
}

#Function that returns legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
###############################
# make Supplementary Figure 6a

why_nrmse<-read.table("../data/SupplementaryFigure6A.txt")
p1<-ggplot(data=why_nrmse, aes(x=why_nrmse_small, y=why_nrmse_big)) + geom_point() + xlab("Squared Difference For RSEM Simulations") + ylab("Squared Difference For Splatter and Polyester Simulations") + xlim(c(0,130)) + ylim(c(0,130)) + theme(text=element_text(size=16), plot.margin = unit(c(3,3,3,3), "lines"))

##############################################
# make supplementary figure 6b

results<-read.table("../data/SupplementaryFigure6B.txt")
results$Simulation<-factor(results$Simulation, levels = c("RSEMsim","bias","unbias"))

results_all<-results[results$ID=="all",]
results_avoid<-results[results$ID=="avoid",]

#make legend
spearmans_graph<-ggplot(data=results_all, aes(x=Tool, y=Value, group=Simulation, fill=factor(Simulation)))  + geom_bar(position="dodge", stat = "summary", fun.y = "mean")
spearmans_graph<- spearmans_graph + scale_x_discrete(labels=c("RSEM", "Salmon Alignment Mode", "Salmon Quasi Mode", "Salmon SMEM Mode", "Sailfish", "eXpress", "Kallisto"))
spearmans_graph<-spearmans_graph + labs(fill="Simulation Method") + scale_fill_manual(labels =c( "RSEM", "Polyester with 3' bias", "Polyester without 3' bias"), values=c("#E69F00", "#56B4E9", "#009E73")) + theme(text=element_text(size=10))
leg<-g_legend(spearmans_graph)

#make graphs
nrmse_graph<-plot_data(results_all, "NRMSE with outliers", "NRMSE", "", c(0,70))
alt_nrmse_graph<-plot_data(results_avoid, "NRMSE with outliers removed", "", "", c(0,70))

ggarrange(p1,ggarrange(nrmse_graph, alt_nrmse_graph, leg, ncol=3, widths=c(2,2,1)),nrow=2, labels=c("A", "B"), heights=c(2,1))
ggsave("../pdfs/SupplementaryFigure6.pdf", plot=last_plot(), height = 6.03 * 1.6, width = 8.66*1.6)
ggsave("../pngs/SupplementaryFigure6.png", plot=last_plot(), height = 6.03 * 1.6, width = 8.66*1.6)


