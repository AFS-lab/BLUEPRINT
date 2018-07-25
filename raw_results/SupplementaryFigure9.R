library(data.table)
library(tidyverse)
##############################################################
# QC

library(ggplot2)
library(ggpubr)
library(reshape2)
library(hydroGOF)
library(DESeq2)
library(SingleCellExperiment)
library(scater)

###############################################################################
# QC FUNCTIONS

MT_reads<-function(counts_data, plates_data){
  #ERCC and MT filtering
  MaleB<-read.table(counts_data, header=T,row.names = 1)
  MaleB<-MaleB[-grep('rRNAfiltered', names(MaleB))]
  
  if (plates_data==TRUE){
    #Load plate data
    plate1<-read.table("../plate1.txt", row.names=1)
    plate2<-read.table("../plate2.txt", row.names = 1)
    
    #Function that formats plate data 
    make_plates<-function(counts_data){
      plates<-colnames(counts_data) %in% rownames(plate1)
      for (i in 1:length(plates)){
        if (plates[i]==TRUE){
          plates[i]<-'batch1'
        }
        else{
          plates[i]<-'batch2'
        }
      }
      return(plates)
    }
    
    batch<-make_plates(MaleB)
  }
  else{
    MaleB<-MaleB[,colnames(MaleB) %in% filter_RSEM$Filename]
    batch<-rep("Batch1", ncol(MaleB))
  }
  
  
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
    feature_controls = list(MT = isSpike(BLUEPRINT_Bs_scater, "MT"))
  )
  
  return(BLUEPRINT_Bs_scater_QC)
}
###############################################################################
# CREATE FILTERS BASED ON QC STATS

#QC script
#Read in QC statistics files
QC_raw<-read.csv("data/raw/read_alignment_qc.csv", header=F)
names(QC_raw)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")

BLUEPRINT_QC_raw<-MT_reads("data/clean_Kallisto_real_Counts.txt", TRUE)
BLUEPRINT_QC_raw<-BLUEPRINT_QC_raw[,BLUEPRINT_QC_raw$pct_counts_MT<10]

#Filters applied based on plots
filter_Kallisto<-QC_raw[QC_raw$NonUnique<350000 & QC_raw$NumReads<4000000 & QC_raw$NumAlignments< 8200000 & QC_raw$Unique<8000000,]
filter_Kallisto<-filter_Kallisto[filter_Kallisto$Filename %in% BLUEPRINT_QC_raw$ids,]
filter_RSEM<-filter_Kallisto

BLUEPRINT_QC_RSEM_sim<-MT_reads("data/clean_ground_truth_counts.txt",FALSE)
BLUEPRINT_QC_RSEM_sim<-BLUEPRINT_QC_RSEM_sim[,BLUEPRINT_QC_RSEM_sim$pct_counts_MT<10]
filter_RSEM<-filter_RSEM[filter_RSEM$Filename %in% BLUEPRINT_QC_RSEM_sim$ids,]

#Read in QC statistics files
QC_sim_poly<-read.csv("../Polyester_Simulations/raw_results/data/simulated/read_alignment_qc.csv", header=F)
names(QC_sim_poly)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")


filter_bias<-QC_sim_poly[grep('bias', QC_sim_poly$Filename),]
filter_bias<-filter_bias[filter_bias$NonUnique<250000,]

filter_unbias<-QC_sim_poly[grep('sample[0-9]', QC_sim_poly$Filename),]
filter_unbias<-filter_unbias[filter_unbias$NonUnique<250000,]

#remove unneeded objects
rm(list=setdiff(ls(), c("filter_RSEM", "filter_bias", "filter_unbias", "filter_Kallisto")))

#########################################################################################
# Format data for ggplot and plot results

get_data_ggplot_ready<-function(path_to_data, colname,filter){
  RSEM_sim<-read.table(path_to_data, header=T,row.names = 1)
  RSEM_sim<-RSEM_sim[,colnames(RSEM_sim) %in% filter$Filename]
  setDT(RSEM_sim, keep.rownames = TRUE)[]
  RSEM_sim<-melt(RSEM_sim)
  means<-group_by(RSEM_sim,rn) %>% summarise(mean(value))
  colnames(means)<-c("rn", colname)
  #means<-means[means$`mean(value)`>0,]
  return(means)
}

Kallisto<-get_data_ggplot_ready("data/clean_Kallisto_real_Counts.txt", "Kallisto",filter_Kallisto)
Kallisto<-Kallisto[order(Kallisto$Kallisto),]
RSEM<-get_data_ggplot_ready("data/clean_ground_truth_counts.txt", "RSEM", filter_RSEM)
RSEM<-RSEM[order(RSEM$RSEM),]
#Splatter<-get_data_ggplot_ready("data/BLUEPRINT_results_benchmark_2/polyester/counts/clean_simulated_counts.txt", "Splatter")

#Create splatter data (have to add zeros)
splatter_sim<-read.table("../Polyester_Simulations/raw_results/data/clean_ground_truth_counts.txt")
RSEM_for_rows<-read.table("data/clean_Kallisto_real_Counts.txt", header=T,row.names = 1)
RSEM_for_rows<-RSEM_for_rows[,1:ncol(splatter_sim)]
colnames(RSEM_for_rows)<-colnames(splatter_sim)
extra_rows<-RSEM_for_rows[!rownames(RSEM_for_rows) %in% rownames(splatter_sim),]
extra_rows[extra_rows!=0]<-0
splatter_sim<-rbind(splatter_sim,extra_rows)

filter_bias$Filename<-sub('samplebias','Cell', filter_bias$Filename)
bias_splatter_sim<-splatter_sim[, colnames(splatter_sim) %in% filter_bias$Filename]
filter_unbias$Filename<-sub('sample','Cell', filter_unbias$Filename)
unbias_splatter_sim<-splatter_sim[, colnames(splatter_sim) %in% filter_unbias$Filename]

rm(splatter_sim)

make_ggplot_results_Splatter<-function(splatter_sim, column){
  setDT(splatter_sim, keep.rownames = TRUE)[]
  splatter_sim<-melt(splatter_sim)
  Splatter<-group_by(splatter_sim,rn) %>% summarise(mean(value))
  colnames(Splatter)<-c("rn", column)
  Splatter<-Splatter[order(Splatter$rn),]
}

bias_splatter_sim<-make_ggplot_results_Splatter(bias_splatter_sim, "Splatter_bias")
unbias_splatter_sim<-make_ggplot_results_Splatter(unbias_splatter_sim, "Splatter_unbias")

summary_df<-merge(Kallisto, RSEM, by="rn")
summary_df<-merge(summary_df, bias_splatter_sim, by="rn")
summary_df<-merge(summary_df, unbias_splatter_sim, by="rn")

write.table(summary_df, gzfile("../figures/data/SupplementaryFigure9.gz"))