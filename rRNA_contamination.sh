#!/bin/bash

filename=$1
path_to_java=$2

./Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir Simulation/indices/STAR --readFilesIn ES_cell_data/$filename"_1.fastq" ES_cell_data/$filename"_2.fastq" --outFileNamePrefix ES_cell_data/$filename --outSAMtype BAM SortedByCoordinate

#Use bam_stat from the RSeQC package to find alignment statistics
source Simulation/venv/bin/activate
split_bam.py -i ES_cell_data/$filename"Aligned.sortedByCoord.out.bam" -r no_chr_mm10_rRNA.bed -o ES_cell_data/$filename
deactivate
rm ES_cell_data/$filename"_1.fastq"
rm ES_cell_data/$filename"_2.fastq"
$path_to_java -XX:MaxHeapSize=1000m -jar ES_cell_data/software/picard.jar SamToFastq \
I=ES_cell_data/$filename".ex.bam" \
FASTQ=ES_cell_data/'rRNA_filtered_'$filename"_1.fastq" \
SECOND_END_FASTQ=ES_cell_data/'rRNA_filtered_'$filename"_2.fastq"
