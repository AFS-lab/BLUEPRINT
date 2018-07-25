library(ggplot2)
library(ggpubr)
library(tidyverse)

##############
# FUNCTIONS

make_ggplot<-function(df, id){
  df<-df[df$ID==id,]
  graph<-ggplot(data=df, aes(x=estimates)) + geom_histogram(binwidth = 10)  +coord_cartesian(ylim=c(0, 100000), xlim=c(0,1500)) + xlab("100 x log2(Ground Truth Expression in TPM + 1)") + ylab("Frequency of Isoforms") 
}

make_zeros_plot<-function(df){
  graph<-ggplot(data=df, aes(y=percentage, x=ID)) + geom_jitter(alpha=0.5, position=position_jitter(width = .2)) + scale_x_discrete(labels=c("Polyester (3' bias)", "Polyester (no bias)", "RSEM"))
  graph<-graph + ylab("% of Isoforms Which Are Unexpressed") + xlab("Simulation") 
}

#########################

ggplot_results<-read.table(gzfile("../data/SupplementaryFigure7.gz"))

p1<-make_ggplot(ggplot_results, "RSEM")
p2<-make_ggplot(ggplot_results, "Polyester_bias")
p3<-make_ggplot(ggplot_results, "Polyester_unbias")

zeros<-ggplot_results %>% group_by(cell, ID) %>% summarise((sum(estimates==0)/128809)*100)
names(zeros)[3]<-c("percentage")
p4<-make_zeros_plot(zeros)

ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, labels=c("A", "B", "C", "D"))
ggsave("../pdfs/SupplementaryFigure7.pdf", plot=last_plot())
ggsave("../pngs/SupplementaryFigure7.png", plot=last_plot())
