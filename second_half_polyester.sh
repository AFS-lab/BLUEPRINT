#!/bin/bash

filename=`echo $1 |awk -F_ '{print $1}'`

./quality_control.sh qualitycontrol $filename Simulation/data/simulated '_1.fq' "simulated"
./quantify.sh Kallisto $filename
./quantify.sh eXpress $filename
./quantify.sh Salmon $filename
./quantify.sh RSEM $filename
./quantify.sh Sailfish $filename
