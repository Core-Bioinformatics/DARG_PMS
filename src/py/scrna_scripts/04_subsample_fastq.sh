#!/bin/bash


new_sample_dir=
sample_name="105_2g_DMSO"
sub_factor=0.4
seed=2025

# go through the fastq files with sample nimse

for i in $(ls ${new_sample_dir}/${sample_name}_*.fastq.gz)
do
    base_name=$(basename $i)
    base_name=${base_name%.fastq.gz}
    first_part=${base_name%_S*_L*_[RI][123]_001}
    second_part=${base_name#$first_part}

    output_path=${new_sample_dir}/${first_part}_sub_${sub_factor}${second_part}
    output_path=${output_path//./_}
    output_path=${output_path}.fastq
    echo $i
    echo $output_path

    seqtk sample -s $seed $i $sub_factor > $output_path
    pigz -p 24 $output_path 
done
