library('splatter')
library('scran')
library('scDD')
library('pscl')

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

###############################################################################
# FUNCTIONS USED BY get_statistics_for_ggplot

#Function which loads data, subsets it and orders it
data_processing<-function(path, filter){
  results<-read.table(path)
  results<-results[,colnames(results) %in% filter$Filename]
  colnames(results)<-sub("bias","", colnames(results))
  colnames(results)<-sub("sample", "Cell", colnames(results))
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  return(results)
}

#Change so load Salmon SMEM counts
#Load Kallisto estimated counts from real data
Kallisto_real_data<-data_processing("data/clean_Kallisto_real_Counts.txt", filter_RSEM)

#Load RSEM simulation counts data
RSEM<-data_processing("data/clean_ground_truth_counts.txt", filter_RSEM)

#Load plate data
plate1<-read.table("../plate1.txt", row.names=1)
plate2<-read.table("../plate2.txt", row.names = 1)

load("../Bs_zeros_removed_25_09_param")

#Function that formats plate data 
make_plates<-function(counts_data){
  plates<-colnames(counts_data) %in% rownames(plate1)
  for (i in 1:length(plates)){
    if (plates[i]==TRUE){
      plates[i]<-'positive'
    }
    else{
      plates[i]<-'negative'
    }
  }
  return(plates)
}

#Function that formats conditions data
make_conditions<-function(counts_data){
  conditions<-colnames(counts_data) %in% colnames(MaleT)
  for (i in 1:length(conditions)){
    if (conditions[i]==TRUE){
      conditions[i]<-1
    } else {
      conditions[i]<-2
    }
  }
  return(conditions)
}

#Function that makes SCEset objects
make_SCE<-function(counts_data){
  #Make phenodata
  Cell<-colnames(counts_data)
  Plate<-make_plates(counts_data)
  pd<-data.frame(Cell,Plate)
  
  #Then make SCESet Object
  rownames(pd) <- pd$Cell
  BLUEPRINT_Bs_scater <- SingleCellExperiment(
    assays = list(counts = as.matrix(counts_data)), 
    colData = pd
  )
  
  return(BLUEPRINT_Bs_scater)
}

simpleParams<-simpleEstimate(data.matrix(Kallisto_real_data, rownames.force=TRUE))
simplesim<-simpleSimulate(simpleParams)

#Perform Lun simulation
lunParams<-lunEstimate(data.matrix(Kallisto_real_data, rownames.force=TRUE))
lunsim<-lunSimulate(lunParams)

lun2sim<-lun2Simulate(lun2Params)

save(simplesim, file="../figures/data/SupplementaryFigure5_simplesim.Rdata")
save(lunsim, file="../figures/data/SupplementaryFigure5_lunsim.Rdata")
save(lun2sim, file="../figures/data/SupplementaryFigure5_lun2sim.Rdata")

RSEM<-make_SCE(RSEM)
save(RSEM, file="../figures/data/SupplementaryFigure5_RSEM.Rdata")

Kallisto<-make_SCE(Kallisto_real_data)
save(Kallisto, file="../figures/data/SupplementaryFigure5_Kallisto.Rdata")

