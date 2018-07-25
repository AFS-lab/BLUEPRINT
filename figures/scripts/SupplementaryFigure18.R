library(ggplot2)
library(ggpubr)
library(reshape2)

##############################
# FUNCTIONS

get_percentage_hist<-function(df){
  no_zeros<-df#[df$value!=0,]
  no_zeros$value<-no_zeros$value==2
  two_iso_percentage<-no_zeros %>% group_by(genes) %>% summarise(percentage(value))
  ggplot(data = two_iso_percentage, aes(x=`percentage(value)`)) + geom_histogram() + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")
} 

#####################################
# MAKE FIGURE S17A

barplot_data<-read.table("../data/Figure6_number_of_isoforms.txt")
B_2_barplot_data<-barplot_data[barplot_data$ID=="B_2",]
B_3_barplot_data<-barplot_data[barplot_data$ID=="B_3",]


B_2_0_1_2_barplot<-ggplot(data=B_2_barplot_data, aes(x=value)) + geom_bar() + xlab("Number of Isoforms Expressed") + ggtitle("Replicate 2")
B_3_0_1_2_barplot<-ggplot(data=B_3_barplot_data, aes(x=value)) + geom_bar() + xlab("Number of Isoforms Expressed") + ggtitle("Replicate 3")

##############################################################
#MAKE FIGURE S17B

exprs_percent_data<-read.table("../data/Figure_6_percent_exprs.txt")
colnames(exprs_percent_data)[2]<-"expression"
colnames(exprs_percent_data)[3]<-"percentage"
B_2_exprs_percent_data<-exprs_percent_data[exprs_percent_data$ID=="B_2",]
B_3_exprs_percent_data<-exprs_percent_data[exprs_percent_data$ID=="B_3",]

B_2_percent_hist<-ggplot(data=B_2_exprs_percent_data, aes(x=percentage)) + geom_histogram() + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")
B_3_percent_hist<-ggplot(data=B_3_exprs_percent_data, aes(x=percentage)) + geom_histogram() + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")

#################################################################
#MAKE FIGURE S17C

B_2_exprs_percent_plot<-ggplot(data =B_2_exprs_percent_data, aes(y=expression,x=percentage)) + geom_point() + ylab("log2(Counts + 1)") + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")
B_3_exprs_percent_plot<-ggplot(data = B_3_exprs_percent_data, aes(y=expression,x=percentage)) + geom_point() + ylab("log2(Counts + 1)") + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")

ggarrange(B_2_0_1_2_barplot, B_3_0_1_2_barplot, B_2_percent_hist, B_3_percent_hist, B_2_exprs_percent_plot, B_3_exprs_percent_plot, ncol =2, nrow=3, labels=c("A", "", "B", "", "C"))
ggsave("../pdfs/SupplementaryFigure18.pdf", plot=last_plot(), height= 225, width=170, units=c("mm"))
ggsave("../pngs/SupplementaryFigure18.png", plot=last_plot(), height= 225, width=170, units=c("mm"))