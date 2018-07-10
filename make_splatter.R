library('splatter')
library('scran')
library('scDD')
library('pscl')
library('limSolve')

###################
#PERFORM QC

single_cell<-read.table("raw_results/data/clean_Kallisto_real_Counts.txt")
single_cell<-single_cell[,colnames(single_cell)!='ERR522956', drop=FALSE]

plate1<-read.table("plate1.txt", row.names=1)
plate2<-read.table("plate2.txt", row.names = 1)

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


ids<-names(single_cell)
batch<-make_plates(single_cell)

anno<-as.data.frame(cbind(batch,ids))
rownames(anno)<-anno$ids

teich_scater <- SingleCellExperiment(
  assays = list(counts = as.matrix(single_cell)), 
  colData = anno
)

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
isSpike(teich_scater, "MT") <- rownames(teich_scater) %in% mt_isoforms

teich_scater_QC <- calculateQCMetrics(
  teich_scater,
  feature_controls = list(MT = isSpike(teich_scater, "MT"))
)

mt_reads<-plotPhenoData(
  teich_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_MT",
             colour = "batch")
)

teich_scater_QC<-teich_scater_QC[,teich_scater_QC$pct_counts_MT<10]

#Read in QC statistics files
QC_raw<-read.csv("raw_results/data/raw/read_alignment_qc.csv", header=FALSE)
names(QC_raw)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")

filter_RSEM<-QC_raw[QC_raw$NonUnique<350000 & QC_raw$NumReads<4000000 & QC_raw$NumAlignments< 8200000 & QC_raw$Unique<8000000,]
filter_RSEM<-filter_RSEM[filter_RSEM$Filename %in% teich_scater_QC$ids,]

rm(list=setdiff(ls(), c("filter_RSEM")))
#######################################
#PERFORM PARAMETER ESTIMATION

#Load Kallisto estimated counts from real data. 
MaleB<-read.table("raw_results/data/clean_Kallisto_real_Counts.txt", header=T,row.names = 1)

#Apply QC filter 
MaleB<-MaleB[,colnames(MaleB) %in% filter_RSEM$Filename]


#Load plate data
plate1<-read.table("plate1.txt", row.names=1)
plate2<-read.table("plate2.txt", row.names = 1)

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

#Only MaleB cells with zeros removed
MaleB<-MaleB[apply(MaleB[,-1], 1, function(x) !all(x==0)),]

#Perform Lun2 simulation parameter estimation - I save this data because the program takes days to run
lun2Params<-lun2Estimate(data.matrix(MaleB, rownames.force=TRUE), as.factor(make_plates(MaleB)), min.size=20)

save(lun2Params, file="Bs_zeros_removed_25_09_param") 

#Use lun2params to simulate gene counts
simulated_counts<-lun2Simulate(params = lun2Params, nGenes=59548)
simulated_counts<-counts(simulated_counts)
write.table(simulated_counts,"raw_results/data/orig_splatter_sim_counts.txt")
