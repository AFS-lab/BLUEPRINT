###############################################################################
# CREATE FILTERS BASED ON QC STATS

#QC script
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(gtable)
library(ggpubr)
library(reshape2)
library(tidyverse)
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
rm(list=setdiff(ls(), c("filter_RSEM")))

################################################################################
#LOAD BLUEPRINT DATA AND PROCESS IT

#Function which loads data, subsets it and orders it
data_processing<-function(path, filter){
  results<-read.table(path)
  results<-results[,colnames(results) %in% filter$Filename]
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  return(results)
}


#load data
ground_truth_B<-data_processing("data/clean_ground_truth_TPM.txt", filter_RSEM)
RSEM_B<-data_processing("data/clean_RSEM_TPM.txt", filter_RSEM)
Salmon_align_B<-data_processing("data/clean_Salmon_align_TPM.txt", filter_RSEM)
Salmon_quasi_B<-data_processing("data/clean_Salmon_quasi_TPM.txt", filter_RSEM)
Salmon_SMEM_B<-data_processing("data/clean_Salmon_SMEM_TPM.txt", filter_RSEM)
Sailfish_B<-data_processing("data/clean_Sailfish_TPM.txt", filter_RSEM)
eXpress_B<-data_processing("data/clean_eXpress_TPM.txt", filter_RSEM)
Kallisto_B<-data_processing("data/clean_Kallisto_TPM.txt", filter_RSEM)


#######################################################################################################

#Remove rows from ground truth matrix with more than 100%, 80%, 60%, 40%, 20% zeros
library(hydroGOF)

#Function to find correlation
correlation_mean<-function(x,y) {
  mean(diag(cor(y,x,method="spearman")))
}

#Function to find standard error in spearmans
correlation_std<-function(x,y) {
  std(diag(cor(y,x,method="spearman")))
}

#Function to find correlation
nrmse_mean<-function(estimates,truth) {
  mean(nrmse(log2(estimates+1),log2(truth+1),))
}

#Function to find standard error
std <- function(x) sd(x)/sqrt(length(x))

#Function to find standard error in spearmans
nrmse_std<-function(estimates,truth) {
  std(nrmse(log2(estimates+1),log2(truth+1)))
}

#Function to remove rows with more than number zeros
remove_zeros<-function(number,truth){
  return(truth[rowSums(truth==0)<(ncol(truth)*number),,drop=FALSE])
}

#Function to keep only rows in ground truth in expression estimates
filter<-function(estimate,truth) {
  return(estimate[rownames(estimate) %in% rownames(truth), , drop=FALSE])
}

#Function that returns misclassification rate
make_misclassify<-function(ground_truth, tool_estimates, threshold_unexpr){
  FN<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates<=threshold_unexpr])
  #print(TP)
  FP<-length(ground_truth[ground_truth<=threshold_unexpr & tool_estimates>threshold_unexpr])
  return((FN+FP)/length(ground_truth))
}

#Function that returns mean misclassification rate per cell
return_mean_misclassify_per_cell<-function(truth_input_data, estimate_input_data){
  
  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_misclassify(truth_input_data[,j], estimate_input_data[,j], 0)
  }
  
  return(mean(do.call(rbind,results)[,1]))
  
}

#Function that returns standard error of misclassification rate per cell
return_std_misclassify_per_cell<-function(truth_input_data, estimate_input_data){
  
  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_misclassify(truth_input_data[,j], estimate_input_data[,j], 0)
  }
  
  return(std(do.call(rbind,results)[,1]))
  
}

tools<-c("RSEM","Salmon Alignment", "Salmon Quasi", "Salmon SMEM", "Sailfish", "eXpress", "Kallisto")

#Write a function which takes a number, ground truth and expression estimates as input, removes rows with more than that percentage of zeros from ground truth, then returns stats
spearmans_rho<-function(number, cell_type, ground_truth,RSEM,Salmon_align, Salmon_quasi, Salmon_SMEM, Sailfish, eXpress, Kallisto){
  #Remove rows with more than number zeros
  ground_truth<-remove_zeros(number,ground_truth)
  
  #Keep only rows in ground truth in expression estimates
  RSEM<-filter(RSEM,ground_truth)
  Salmon_align<-filter(Salmon_align,ground_truth)
  Salmon_quasi<-filter(Salmon_quasi,ground_truth)
  Salmon_SMEM<-filter(Salmon_SMEM,ground_truth)
  Sailfish<-filter(Sailfish,ground_truth)
  eXpress<-filter(eXpress,ground_truth)
  Kallisto<-filter(Kallisto,ground_truth)
  
  #create empty vectors for results
  spearmans_results<- vector(mode="numeric", length=0)
  spearmans_error<- vector(mode="numeric", length=0)
  nrmse_results<-vector(mode="numeric", length=0)
  nrmse_error<-vector(mode="numeric", length=0)
  misclassify_results<-vector(mode="numeric", length=0)
  misclassify_error<-vector(mode="numeric", length=0)
  
  
  #Find spearmans rho
  spearmans_results[1]<-correlation_mean(RSEM,ground_truth)
  spearmans_results[2]<-correlation_mean(Salmon_align,ground_truth)
  spearmans_results[3]<-correlation_mean(Salmon_quasi,ground_truth)
  spearmans_results[4]<-correlation_mean(Salmon_SMEM, ground_truth)
  spearmans_results[5]<-correlation_mean(Sailfish,ground_truth)
  spearmans_results[6]<-correlation_mean(eXpress,ground_truth)
  spearmans_results[7]<-correlation_mean(Kallisto,ground_truth)
  
  #find standard error
  spearmans_error[1]<-correlation_std(RSEM,ground_truth)
  spearmans_error[2]<-correlation_std(Salmon_align,ground_truth)
  spearmans_error[3]<-correlation_std(Salmon_quasi,ground_truth)
  spearmans_error[4]<-correlation_std(Salmon_SMEM, ground_truth)
  spearmans_error[5]<-correlation_std(Sailfish,ground_truth)
  spearmans_error[6]<-correlation_std(eXpress,ground_truth)
  spearmans_error[7]<-correlation_std(Kallisto,ground_truth)
  
  #Find nrmse
  nrmse_results[1]<-nrmse_mean(RSEM,ground_truth)
  nrmse_results[2]<-nrmse_mean(Salmon_align,ground_truth)
  nrmse_results[3]<-nrmse_mean(Salmon_quasi,ground_truth)
  nrmse_results[4]<-nrmse_mean(Salmon_SMEM, ground_truth)
  nrmse_results[5]<-nrmse_mean(Sailfish,ground_truth)
  nrmse_results[6]<-nrmse_mean(eXpress,ground_truth)
  nrmse_results[7]<-nrmse_mean(Kallisto,ground_truth)
  
  #return vector of spearmans rho values
  nrmse_error[1]<-nrmse_std(RSEM,ground_truth)
  nrmse_error[2]<-nrmse_std(Salmon_align,ground_truth)
  nrmse_error[3]<-nrmse_std(Salmon_quasi,ground_truth)
  nrmse_error[4]<-nrmse_std(Salmon_SMEM, ground_truth)
  nrmse_error[5]<-nrmse_std(Sailfish,ground_truth)
  nrmse_error[6]<-nrmse_std(eXpress,ground_truth)
  nrmse_error[7]<-nrmse_std(Kallisto,ground_truth)
  
  #Find nrmse
  #misclassify_results[1]<-return_mean_misclassify_per_cell(RSEM,ground_truth)
  #misclassify_results[2]<-return_mean_misclassify_per_cell(Salmon_align,ground_truth)
  #misclassify_results[3]<-return_mean_misclassify_per_cell(Salmon_quasi,ground_truth)
  #misclassify_results[4]<-return_mean_misclassify_per_cell(Salmon_SMEM, ground_truth)
  #misclassify_results[5]<-return_mean_misclassify_per_cell(Sailfish,ground_truth)
  #misclassify_results[6]<-return_mean_misclassify_per_cell(eXpress,ground_truth)
  #misclassify_results[7]<-return_mean_misclassify_per_cell(Kallisto,ground_truth)
  
  #misclassify_error[1]<-return_std_misclassify_per_cell(RSEM,ground_truth)
  #misclassify_error[2]<-return_std_misclassify_per_cell(Salmon_align,ground_truth)
  #misclassify_error[3]<-return_std_misclassify_per_cell(Salmon_quasi,ground_truth)
  #misclassify_error[4]<-return_std_misclassify_per_cell(Salmon_SMEM, ground_truth)
  #misclassify_error[5]<-return_std_misclassify_per_cell(Sailfish,ground_truth)
  #misclassify_error[6]<-return_std_misclassify_per_cell(eXpress,ground_truth)
  #misclassify_error[7]<-return_std_misclassify_per_cell(Kallisto,ground_truth)
  
  return(cbind(tools,spearmans_results, spearmans_error, nrmse_results, nrmse_error, rep(number, 7), rep(cell_type, 7)))
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

#make suitable dataframe for ggplot
hundred_B<-spearmans_rho(1, "B", ground_truth_B, RSEM_B,Salmon_align_B, Salmon_quasi_B, Salmon_SMEM_B, Sailfish_B, eXpress_B, Kallisto_B)
eighty_B<-spearmans_rho(0.8,"B",ground_truth_B, RSEM_B,Salmon_align_B, Salmon_quasi_B, Salmon_SMEM_B, Sailfish_B, eXpress_B, Kallisto_B)
sixty_B<-spearmans_rho(0.6,"B",ground_truth_B, RSEM_B,Salmon_align_B, Salmon_quasi_B, Salmon_SMEM_B, Sailfish_B, eXpress_B, Kallisto_B)
forty_B<-spearmans_rho(0.4,"B",ground_truth_B, RSEM_B,Salmon_align_B, Salmon_quasi_B, Salmon_SMEM_B, Sailfish_B, eXpress_B, Kallisto_B)
twenty_B<-spearmans_rho(0.2,"B",ground_truth_B, RSEM_B,Salmon_align_B, Salmon_quasi_B, Salmon_SMEM_B, Sailfish_B, eXpress_B, Kallisto_B)

figure_5a<-rbind(hundred_B, eighty_B, sixty_B, forty_B,twenty_B)
colnames(figure_5a)[6]<-"percentage_zeros"
figure_5a[,6]<-as.numeric(figure_5a[,6]) * 100
colnames(figure_5a)[7]<-"cell_type"

write.table(figure_5a,"../figures/data/Figure5a.txt")
