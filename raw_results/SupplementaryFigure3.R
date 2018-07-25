#QC script
library(ggplot2)
library(ggpubr)
library(scater, quietly = TRUE)
library(knitr)
library(SingleCellExperiment)
options(stringsAsFactors = FALSE)

QC_sim_poly<-read.csv("../Polyester_Simulations/raw_results/data/simulated/read_alignment_qc.csv", header=F)
names(QC_sim_poly)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")

#Extract biased cells
QC_bias<-QC_sim_poly[grep('bias', QC_sim_poly$Filename),]

write.table(QC_bias, "../figures/data/SupplementaryFigure3_reads_alignment_data.txt")

MaleB<-read.table("../Polyester_Simulations/raw_results/data/clean_ground_truth_counts.txt", header=T,row.names = 1)

batch<-rep("Batch1", ncol(MaleB))
ids<-colnames(MaleB)

anno<-as.data.frame(cbind(batch,ids))
rownames(anno)<-anno$ids
BLUEPRINT_Bs_scater <- SingleCellExperiment(
  assays = list(counts = as.matrix(MaleB)), 
  colData = anno
)

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
isSpike(BLUEPRINT_Bs_scater, "MT") <- rownames(BLUEPRINT_Bs_scater) %in% mt_isoforms


BLUEPRINT_Bs_scater_QC <- calculateQCMetrics(
  BLUEPRINT_Bs_scater,
  feature_controls = list(MT = isSpike(BLUEPRINT_Bs_scater, "MT")))


save(BLUEPRINT_Bs_scater_QC, file="../figures/data/SupplementaryFigure3_scater_object.RData")
