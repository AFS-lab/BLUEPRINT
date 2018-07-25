library(ggplot2)
library(ggpubr)
library(reshape2)
library(hydroGOF)
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
rm(list=setdiff(ls(), c("filter_RSEM", "filter_bias", "filter_unbias")))

#For now I won't apply any filters on the simulated data

###############################################################################
# FUNCTIONS USED BY get_statistics_for_ggplot

#Function which loads data, subsets it and orders it
data_processing<-function(path, filter){
  results<-read.table(path)
  results<-results[,colnames(results) %in% filter$Filename]
  colnames(results)<-sub("sample", "Cell", colnames(results))
  colnames(results)<-sub("bias","", colnames(results))
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  return(results)
}

#Function to find correlation
correlation<-function(x,y) {
  (diag(cor(y,x,method="spearman")))
}

#Function that returns precision
make_precision<-function(ground_truth, tool_estimates, threshold_unexpr){
  TP<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates>threshold_unexpr])
  #print(TP)
  FP<-length(ground_truth[ground_truth<=threshold_unexpr & tool_estimates>threshold_unexpr])
  return(TP/(TP+FP))
}

#Function that returns precision value per cell
return_precision_per_cell<-function(truth_input_data, estimate_input_data){
  
  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_precision(truth_input_data[,j], estimate_input_data[,j], 0)
  }
  
  return(do.call(rbind,results))
  
}


#Function that returns recall
make_recall<-function(ground_truth, tool_estimates, threshold_unexpr){
  TP<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates>threshold_unexpr])
  #print(TP)
  FN<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates<=threshold_unexpr])
  return(TP/(TP+FN))
}

#Function that returns recall per cell
return_recall_per_cell<-function(truth_input_data, estimate_input_data){
  
  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_recall(truth_input_data[,j], estimate_input_data[,j], 0)
  }
  
  return(do.call(rbind,results))
  
}

#Function that returns F1
find_F1<-function(precision,recall){
  F1<-2*((precision*recall)/(precision + recall))
  return(F1)
}

#Function that returns legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

#Function which returns name of transcripts with (s-o)^2>30 for splatter/polyester 
shrink_nrmse<-function(){
  RSEM<-data_processing("../Polyester_Simulations/raw_results/data/clean_RSEM_TPM.txt", filter_bias)
  ground_truth<-read.table("../Polyester_Simulations/raw_results/data/clean_ground_truth_TPM.txt")
  ground_truth<-ground_truth[,colnames(ground_truth) %in% colnames(RSEM)]
  extra_rows<-RSEM[!rownames(RSEM) %in% rownames(ground_truth),]
  extra_rows[extra_rows!=0]<-0
  ground_truth<-rbind(ground_truth,extra_rows)
  ground_truth<-ground_truth[ , order(colnames(ground_truth))]
  ground_truth<-ground_truth[order(rownames(ground_truth)),]
  why_nrmse_big<-(log2(RSEM+1)-log2(ground_truth+1))^2
  
  why_nrmse_big<-rowMeans(why_nrmse_big)
  
  RSEM_r<-data_processing("data/clean_RSEM_TPM.txt", filter_RSEM)
  ground_truth_r<-data_processing("data/clean_ground_truth_TPM.txt", filter_RSEM)
  why_nrmse_small<-(log2(RSEM_r+1)-log2(ground_truth_r+1))^2
  why_nrmse_small<-rowMeans(why_nrmse_small)
  plot(why_nrmse_small, why_nrmse_big, ylim=c(0,130), xlim=c(0,130))
  problem_isoforms<-why_nrmse_big[why_nrmse_big>30]
  return(names(why_nrmse_big[why_nrmse_big>30]))
}

find_nrmse<-function(ground_truth, tool, rows_to_avoid){
  ground_truth<-ground_truth[!rownames(ground_truth) %in% rows_to_avoid,]
  tool<-tool[!rownames(tool) %in% rows_to_avoid,]
  nrmse<-nrmse(log2(tool+1),log2(ground_truth+1))
  return(nrmse)
}


################################################################################
#LOAD DATA AND PROCESS IT

#big function which loads data, processes it, finds statistics and returns it. Note this function assumes your results matrices are in a single directory and there are no other text files in that directory
get_statistics_for_ggplot<-function(path_to_dir,filter,add_isoforms, nrmse){
  #TODO: For given path, open + process all files in directory
  filenames <- list.files(path_to_dir, pattern="*_TPM.txt", full.names=TRUE)
  
  #create objects for files
  for (i in 1:length(filenames)){
    program<-strsplit(filenames[i],"clean_")[[1]][2]
    program<-strsplit(program,"_TPM.txt")[[1]][1]
    print(program)
    assign(program, data_processing(filenames[i], filter))
    #print(head(get(program)))
  }
  
  #If not RSEM, need to add unexpressed isoforms which were not input into the simulation process
  if (add_isoforms==TRUE){
    ground_truth<-read.table("../Polyester_Simulations/raw_results/data/clean_ground_truth_TPM.txt")
    ground_truth<-ground_truth[,colnames(ground_truth) %in% colnames(Salmon_align)]
    extra_rows<-Salmon_align[!rownames(Salmon_align) %in% rownames(ground_truth),]
    extra_rows[extra_rows!=0]<-0
    ground_truth<-rbind(ground_truth,extra_rows)
    ground_truth<-ground_truth[ , order(colnames(ground_truth))]
    ground_truth<-ground_truth[order(rownames(ground_truth)),]
    #return(ground_truth)
  }
  
  
  if (nrmse==TRUE){
    rows_to_avoid<-shrink_nrmse()
    print(rows_to_avoid)
    
    #Find NRMSE for each method
    RSEM_nmrse<-find_nrmse(ground_truth,RSEM,rows_to_avoid)
    Salmon_align_nmrse<-find_nrmse(ground_truth,Salmon_align,rows_to_avoid)
    Salmon_quasi_nmrse<-find_nrmse(ground_truth,Salmon_quasi,rows_to_avoid)
    Salmon_SMEM_nmrse<-find_nrmse(ground_truth,Salmon_SMEM,rows_to_avoid)
    Sailfish_nmrse<-find_nrmse(ground_truth,Sailfish,rows_to_avoid)
    eXpress_nmrse<-find_nrmse(ground_truth,eXpress,rows_to_avoid)
    Kallisto_nmrse<-find_nrmse(ground_truth,Kallisto,rows_to_avoid)
    
  } else {
    
    RSEM_nmrse<-nrmse(log2(RSEM+1), log2(ground_truth +1))
    Salmon_align_nmrse<-nrmse(log2(Salmon_align+1), log2(ground_truth +1))
    Salmon_quasi_nmrse<-nrmse(log2(Salmon_quasi+1), log2(ground_truth +1))
    Salmon_SMEM_nmrse<-nrmse(log2(Salmon_SMEM+1), log2(ground_truth +1))
    Sailfish_nmrse<-nrmse(log2(Sailfish+1), log2(ground_truth +1))
    eXpress_nmrse<-nrmse(log2(eXpress+1), log2(ground_truth +1))
    Kallisto_nmrse<-nrmse(log2(Kallisto+1), log2(ground_truth +1))
  }
  
  
  #Store in a dataframe
  nrmse_data<-melt(rbind(RSEM_nmrse,Salmon_align_nmrse, Salmon_quasi_nmrse, Salmon_SMEM_nmrse, Sailfish_nmrse, eXpress_nmrse, Kallisto_nmrse))
  #nrmse_data<-melt(rbind(Salmon_align_nmrse, Salmon_quasi_nmrse, Salmon_SMEM_nmrse, Sailfish_nmrse, eXpress_nmrse, Kallisto_nmrse))
  nrmse_data<-cbind("nrmse", nrmse_data)
  nrmse_data<-nrmse_data[,colnames(nrmse_data)!="Var2"]
  colnames(nrmse_data)<-c("Statistic","Tool","Value")
  print(head(nrmse_data))
  
  return(nrmse_data)
}


ggplot_results_RSEM<-get_statistics_for_ggplot("data", filter_RSEM,FALSE, TRUE)
ggplot_results_unbias<-get_statistics_for_ggplot("../Polyester_Simulations/raw_results/data", filter_unbias,TRUE,TRUE)
ggplot_results_bias<-get_statistics_for_ggplot("../Polyester_Simulations/raw_results/data", filter_bias,TRUE,TRUE)

ggplot_results_RSEM<-cbind("RSEMsim", ggplot_results_RSEM)
ggplot_results_unbias<-cbind("unbias", ggplot_results_unbias)
ggplot_results_bias<-cbind("bias", ggplot_results_bias)

colnames(ggplot_results_RSEM)<-c("Simulation", "Statistic", "Tool", "Value")
colnames(ggplot_results_bias)<-colnames(ggplot_results_RSEM)
colnames(ggplot_results_unbias)<-colnames(ggplot_results_RSEM)

results_avoid<-rbind(ggplot_results_RSEM, ggplot_results_bias, ggplot_results_unbias)
results_avoid$Simulation<-factor(results_avoid$Simulation, levels = c("RSEMsim","bias","unbias"))

rm(list=ls(pattern="ggplot_results"))

ggplot_results_RSEM<-get_statistics_for_ggplot("data", filter_RSEM,FALSE, FALSE)
ggplot_results_unbias<-get_statistics_for_ggplot("../Polyester_Simulations/raw_results/data", filter_unbias,TRUE,FALSE)
ggplot_results_bias<-get_statistics_for_ggplot("../Polyester_Simulations/raw_results/data", filter_bias,TRUE,FALSE)

ggplot_results_RSEM<-cbind("RSEMsim", ggplot_results_RSEM)
ggplot_results_unbias<-cbind("unbias", ggplot_results_unbias)
ggplot_results_bias<-cbind("bias", ggplot_results_bias)

colnames(ggplot_results_RSEM)<-c("Simulation", "Statistic", "Tool", "Value")
colnames(ggplot_results_bias)<-colnames(ggplot_results_RSEM)
colnames(ggplot_results_unbias)<-colnames(ggplot_results_RSEM)

results_all<-rbind(ggplot_results_RSEM, ggplot_results_bias, ggplot_results_unbias)
results_all$Simulation<-factor(results_all$Simulation, levels = c("RSEMsim","bias","unbias"))

results_all<-cbind(results_all, ID=rep("all", nrow(results_all)))
results_avoid<-cbind(results_avoid, ID=rep("avoid", nrow(results_avoid)))

results<-rbind(results_all, results_avoid)
write.table(results, "../figures/data/SupplementaryFigure6B.txt")

#Plot explaining why difference in NRMSE 
RSEM<-data_processing("../Polyester_Simulations/raw_results/data/clean_RSEM_TPM.txt", filter_bias)
ground_truth<-read.table("../Polyester_Simulations/raw_results/data/clean_ground_truth_TPM.txt")
ground_truth<-ground_truth[,colnames(ground_truth) %in% colnames(RSEM)]
extra_rows<-RSEM[!rownames(RSEM) %in% rownames(ground_truth),]
extra_rows[extra_rows!=0]<-0
ground_truth<-rbind(ground_truth,extra_rows)
ground_truth<-ground_truth[ , order(colnames(ground_truth))]
ground_truth<-ground_truth[order(rownames(ground_truth)),]
why_nrmse_big<-(log2(RSEM+1)-log2(ground_truth+1))^2
why_nrmse_big<-rowMeans(why_nrmse_big)

RSEM_r<-data_processing("data/clean_RSEM_TPM.txt", filter_RSEM)
ground_truth_r<-data_processing("data/clean_ground_truth_TPM.txt", filter_RSEM)
why_nrmse_small<-(log2(RSEM_r+1)-log2(ground_truth_r+1))^2
why_nrmse_small<-rowMeans(why_nrmse_small)
why_nrmse<-data.frame(why_nrmse_big, why_nrmse_small)

write.table(why_nrmse, "../figures/data/SupplementaryFigure6A.txt")
