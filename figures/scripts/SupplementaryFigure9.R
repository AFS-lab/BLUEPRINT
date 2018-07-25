library(ggplot2)
library(ggpubr)

#load data
summary_df<-read.table(gzfile("../data/SupplementaryFigure9.gz"))

summary_df$rn <- factor(summary_df$rn, levels = summary_df$rn[order(summary_df$Kallisto)])

#plot data
p1<-ggplot(data=summary_df, aes(x=Kallisto, y=Splatter_bias)) + geom_point() + xlab("Kallisto Counts") + ylab("Splatter Counts (3' bias)")
p2<-ggplot(data=summary_df, aes(x=Kallisto, y=Splatter_unbias)) + geom_point() + xlab("Kallisto Counts") + ylab("Splatter Counts (no coverage bias)")
p3<-ggplot(data=summary_df, aes(x=Kallisto, y=RSEM)) + geom_point() + xlab("Kallisto Counts") + ylab("RSEM Counts")

ggarrange(p1,p2,p3, ncol=2, nrow=2, labels=c("A","B", "C"))
ggsave("../pdfs/SupplementaryFigure9.pdf", plot=last_plot())
ggsave("../pngs/SupplementaryFigure9.png")
