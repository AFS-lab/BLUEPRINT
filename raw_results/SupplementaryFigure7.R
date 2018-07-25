library(ggplot2)
library(ggpubr)
library(reshape2)
library(hydroGOF)
library(scater)
library(SingleCellExperiment)
library(data.table)
library(tidyverse)
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
filter_RSEM<-QC_raw[QC_raw$NonUnique<350000 & QC_raw$NumReads<4000000 & QC_raw$NumAlignments< 8200000 & QC_raw$Unique<8000000,]
filter_RSEM<-filter_RSEM[filter_RSEM$Filename %in% BLUEPRINT_QC_raw$ids,]

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
rm(list=setdiff(ls(), c("filter_RSEM", "filter_bias", "filter_unbias")))

################################################################################
#LOAD DATA AND PROCESS IT

#Function which loads data, subsets it and orders it
data_processing<-function(path, filter,RSEM){
  results<-read.table(path)
  if (RSEM==TRUE){
    results<-results[,colnames(results) %in% filter$Filename]
    results<-results[ , order(colnames(results))]
    results<-results[order(rownames(results)),]
    return(results)
  }
  else{
    #fix difference in cell names between filter and table
    filter$Filename<-sub("sample", "Cell", filter$Filename)
    filter$Filename<-sub("bias","", filter$Filename)
    results<-results[,colnames(results) %in% filter$Filename]
    results<-results[ , order(colnames(results))]
    results<-results[order(rownames(results)),]
    
    #add zeros
    RSEM_for_rows<-read.table("data/clean_Kallisto_real_Counts.txt", header=T,row.names = 1)
    RSEM_for_rows<-RSEM_for_rows[,1:ncol(results)]
    colnames(RSEM_for_rows)<-colnames(results)
    extra_rows<-RSEM_for_rows[!rownames(RSEM_for_rows) %in% rownames(results),]
    extra_rows[extra_rows!=0]<-0
    results<-rbind(results,extra_rows)
    
    return(results)
  }
  
}

make_ggplot<-function(path,filter,RSEM, id){
  #load data
  ground_truth<-data_processing(path, filter,RSEM)
  setDT(ground_truth, keep.rownames = TRUE)[]
  
  cells<-colnames(ground_truth)
  
  ground_truth_expr<-ground_truth %>% gather(cells[2:ncol(ground_truth)], key="cell", value="estimates")
  ground_truth_expr$estimates<-100*log2(ground_truth_expr$estimates +1)
  ground_truth_expr<-cbind(ground_truth_expr, ID=rep(id, nrow(ground_truth_expr)))
  return(ground_truth_expr)
}

p1<-make_ggplot("data/clean_ground_truth_TPM.txt", filter_RSEM,TRUE, "RSEM")
p2<-make_ggplot("../Polyester_Simulations/raw_results/data/clean_ground_truth_TPM.txt",filter_bias,FALSE, "Polyester_bias")
p3<-make_ggplot("../Polyester_Simulations/raw_results/data/clean_ground_truth_TPM.txt", filter_unbias,FALSE, "Polyester_unbias")

results<-rbind(p1,p2,p3)

write.table(results,gzfile("../figures/data/SupplementaryFigure7.gz"))
