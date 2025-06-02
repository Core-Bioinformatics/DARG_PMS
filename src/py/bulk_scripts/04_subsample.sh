#!/bin/bash

source_folder=
mkdir -p ${source_folder}/04.Subsampled
cp -r ${source_folder}/03.LanesMerged/* ${source_folder}/04.Subsampled

seqtk sample -s 2025 ${source_folder}/04.Subsampled/C1/C1_L001_R1_001.fastq.gz 20000000 > ${source_folder}/04.Subsampled/C1/C1_L001_R1_001.fastq
pigz -p 24 ${source_folder}/04.Subsampled/C1/C1_L001_R1_001.fastq
seqtk sample -s 2025 ${source_folder}/04.Subsampled/C1/C1_L001_R2_001.fastq.gz 20000000 > ${source_folder}/04.Subsampled/C1/C1_L001_R2_001.fastq
pigz -p 24 ${source_folder}/04.Subsampled/C1/C1_L001_R2_001.fastq
