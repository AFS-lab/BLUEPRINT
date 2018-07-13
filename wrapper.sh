# #!/bin/bash
#
# #rename command line arguments
# path_to_ref_fasta=$1
# path_to_ref_gtf=$2
#
# #Make indexes + RSEM ref for RSEM simulations
# ./setup.sh setup
# ./RSEM_ref.sh make_ref $path_to_ref_gtf $path_to_ref_fasta
# ./make_indexes.sh $path_to_ref_gtf $path_to_ref_fasta
# gunzip ES_cell_data/*
#
# #Make indexes + RSEM ref for Splatter/Polyester Simulations
# cd Polyester_Simulations
# ./RSEM_ref.sh make_ref $path_to_ref_gtf $path_to_ref_fasta
# ./make_indexes.sh $path_to_ref_gtf $path_to_ref_fasta
#
# #Run Kallisto + calculate QC stats for Splatter simulations
# for i in ../ES_cell_data/*_1.fastq;
# do
#   ./first_half_polyester.sh $i
# done
#
# #Format data + move into correct dirs
# python generate.py Kallisto_real `pwd` Simulation/Kallisto_results_real_data
# chmod +x Kallisto_real_Counts.sh
# ./Kallisto_real_Counts.sh
# ./clean_data.sh
# mkdir raw_results/data
# cp -r Simulation/QC_stats/raw raw_results/data/
# cp Simulation/results_matrices/clean* raw_results/data/
#
# #Run splatter + polyester simulations
# Rscript make_splatter.R
# ./control_polyester_script.sh
# ./convert_fasta_to_fastq.sh

#Do benchmarking on splatter/polyester simulated data
for i in Simulation/data/simulated/*_1.fq;
do
  filename=`echo $i | awk -F/ '{print $NF}'`
  bsub -n8 -R"span[hosts=1]" -c 99999 -G team_hemberg -q normal -o $TEAM/temp.logs/output.$filename -e $TEAM/temp.logs/error.$filename -R"select[mem>100000] rusage[mem=100000]" -M100000 ./second_half_polyester.sh $filename
done

#Make results matrices and move to appropriate location
./make_matrix.sh make_matrix RSEM
./make_matrix.sh make_matrix eXpress
./make_matrix.sh make_matrix Kallisto
./make_matrix.sh make_matrix Sailfish
./make_matrix.sh make_matrix Salmon
./make_matrix.sh make_matrix ground_truth
python generate.py Kallisto_real `pwd` Simulation/Kallisto_results_real_data
chmod +x Kallisto_real_Counts.sh
./Kallisto_real_Counts.sh
./clean_data.sh

#move results into appropriate location
cp Simulation/results_matrices/clean* raw_results/data/
cp -r Simulation/QC_stats/raw raw_results/data/
cp -r Simulation/QC_stats/simulated raw_results/data/

#Now the RSEM simulations...
cd ../ES_cell_data

#Download bulk data
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR523/SRR5237781/SRR5237781.sra
./software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --split-3 SRR5237781.sra

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR523/SRR5237782/SRR5237782.sra
./software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --split-3 SRR5237782.sra

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR523/SRR5237783/SRR5237783.sra
./software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --split-3 SRR5237783.sra


#Simulate + benchmark RSEM simulated data
for i in ES_cell_data/*_1.fastq;
do
  num_jobs=`bjobs | wc -l`
  max_jobs=30
  filename=`echo $i | awk -F/ '{print $2}'`

  #This prevents the number of queued jobs greatly exceeding 30.
  while [[ $num_jobs -gt $max_jobs ]];
  do
    sleep 100
    num_jobs=`bjobs | wc -l`
  done

  bsub -n8 -R"span[hosts=1]" -c 99999 -G team_hemberg -q normal -o $TEAM/temp.logs/output.$filename -e $TEAM/temp.logs/error.$filename -R"select[mem>100000] rusage[mem=100000]" -M100000 ./cell_level_analysis.sh $filename
done

#make results matrices
./make_matrix.sh make_matrix RSEM
./make_matrix.sh make_matrix eXpress
./make_matrix.sh make_matrix Kallisto
./make_matrix.sh make_matrix Sailfish
./make_matrix.sh make_matrix Salmon
./make_matrix.sh make_matrix ground_truth
python generate.py Kallisto_real `pwd` Simulation/Kallisto_results_real_data
chmod +x Kallisto_real_Counts.sh
./Kallisto_real_Counts.sh
./clean_data.sh

#move results into appropriate location
cp Simulation/results_matrices/clean* raw_results/data/
cp -r Simulation/QC_stats/raw raw_results/data/
cp -r Simulation/QC_stats/simulated raw_results/data/

#format data to make figures
cd raw_results
Rscript Figure2.R
Rscript Figure4.R
Rscript Figure5a.R
Rscript Figure5b.R
Rscript Figure6.R
Rscript SupplementaryFigure10.R
Rscript SupplementaryFigure11.R
Rscript SupplementaryFigure12.R

#make figure pdfs
cd ../figures/scripts
Rscript Figure2.R
Rscript Figure5.R
Rscript Figure6.R
Rscript SupplementaryFigure10.R
Rscript SupplementaryFigure11.R
Rscript SupplementaryFigure12.R
Rscript SupplementaryFigure16.R
Rscript SupplementaryFigure17.R



python generate.py Kallisto_real `pwd` Simulation/Kallisto_results_real_data
chmod +x Kallisto_real_Counts.sh
./Kallisto_real_Counts.sh
./clean_data.sh

Rscript make_splatter.R
