#!/bin/bash


filename=`echo $1 | awk -F/ '{print $3}' | awk -F_ '{print $1}'`

./quality_control.sh qualitycontrol $filename ../ES_cell_data '_1.fastq' 'raw'
./quantify_real_data.sh Kallisto $filename ../ES_cell_data
