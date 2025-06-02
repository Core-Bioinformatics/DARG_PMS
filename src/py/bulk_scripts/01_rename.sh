#!/bin/bash

source_folder="01.RawData"
target_folder="03.Renamed"

for sample_name in $(ls $source_folder); do
    echo "Renaming $sample_name"
    i=1
    mkdir -p $target_folder/$sample_name
    for fl in $(ls $source_folder/$sample_name/*1.fq.gz | rev | sort | rev); do
        new_name="${sample_name}_L00${i}_R1_001.fastq.gz"
        cp $fl $target_folder/$sample_name/$new_name
        f1=$(echo $fl | sed 's/1.fq/2.fq/')
        new_name="${sample_name}_L00${i}_R2_001.fastq.gz"
        cp $f1 $target_folder/$sample_name/$new_name
        i=$((i+1))
    done
done
