#!/bin/bash

cd Simulation/data

for file in simulated*[0-9];
do
  number=`echo $file| awk -Fd '{print $2}'`
  cp $file/sample_01_1.fasta simulated/sample_$number"_1.fasta"
  cp $file/sample_01_2.fasta simulated/sample_$number"_2.fasta"
done

for file in sim_bias*[0-9];
do
  number=`echo $file| awk -Fd '{print $2}'`
  cp $file/sample_01_1.fasta simulated/samplebias_$number"_1.fasta"
  cp $file/sample_01_2.fasta simulated/samplebias_$number"_2.fasta"
done


cd ../..
for file in Simulation/data/simulated/sample*_1.fasta;
do
  filepath=`echo $file | awk -F_ '{print $1"_"$2}'`
  ./Simulation/bbmap/reformat.sh in1=$filepath"_1.fasta" in2=$filepath"_2.fasta" out1=$filepath"_1.fq" out2=$filepath"_2.fq" qfake=30
done
