#!/bin/bash

filename=`echo $1 |awk -F_ '{print $1}'`

#There is a bug in Polyester which produces 2 reads with the same name for each million and first read. This upsets RSEM so I rename the second of each of these reads here:
sed -i 's/read1000001/read1000001a/2' Simulation/data/simulated/$1
sed -i 's/read1000001/read1000001a/2' Simulation/data/simulated/$filename'_2.fq'
sed -i 's/read2000001/read2000001a/2' Simulation/data/simulated/$1
sed -i 's/read2000001/read2000001a/2' Simulation/data/simulated/$filename'_2.fq'
sed -i 's/read3000001/read3000001a/2' Simulation/data/simulated/$1
sed -i 's/read3000001/read3000001a/2' Simulation/data/simulated/$filename'_2.fq'
sed -i 's/read4000001/read4000001/2' Simulation/data/simulated/$1
sed -i 's/read4000001/read4000001a/2' Simulation/data/simulated/$filename'_2.fq'
sed -i 's/read5000001/read5000001a/2' Simulation/data/simulated/$1
sed -i 's/read5000001/read5000001a/2' Simulation/data/simulated/$filename'_2.fq'
# ./quality_control.sh qualitycontrol $filename Simulation/data/simulated '_1.fq' "simulated"
# ./quantify.sh Kallisto $filename
# ./quantify.sh eXpress $filename
# ./quantify.sh Salmon $filename
./quantify.sh RSEM $filename
# ./quantify.sh Sailfish $filename
