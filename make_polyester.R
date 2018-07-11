library('polyester')

#read and name command line arguments
args <- commandArgs(trailingOnly = TRUE)
number<-as.numeric(args[1])

simulated_counts<-read.table("raw_results/data/simulated_counts.txt")
outdirectory<-paste("Simulation/data/simulated",number, sep="")
outdirectory_bias<-paste("Simulation/data/sim_bias",number, sep="")

#execute polyester simulation with no positional bias
simulate_experiment_countmat("polyester_transcripts", readmat = matrix(simulated_counts[,number]), outdir = outdirectory, paired = TRUE, seed = 0, distr='normal', error_model='illumina4', bias='none')
simulate_experiment_countmat("polyester_transcripts", readmat = matrix(simulated_counts[,number]), outdir = outdirectory_bias, paired = TRUE, seed = 0, distr='normal', error_model='illumina4', bias='cdnaf')
