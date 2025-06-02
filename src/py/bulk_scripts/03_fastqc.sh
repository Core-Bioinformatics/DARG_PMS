#!/bin/bash

source_folder="03.LanesMerged"

fastq_list=$(ls ${source_folder}/*/*R1*)
echo $fastq_list
fastqc -t 30 $fastq_list