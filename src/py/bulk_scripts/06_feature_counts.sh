#!/bin/bash

source_folder=
sourceDir="${source_folder}/05.Star_mapped"
targetDir="${source_folder}/06.Feature_counts"
gtfFile="GRch38_p113.gtf"
nThreads=16

if [ ! -d $targetDir ]; then
    mkdir -p $targetDir
fi

for sample in $(ls $sourceDir); do
    # check if it starts with multiqc
    if [[ $sample == FASTQC_sub_multiqc* ]]; then
        continue
    fi

    for f in $(ls $sourceDir/$sample/*bam); do
        prefix=$(echo $f | rev | cut -d / -f 1 | rev | sed 's/_Aligned.sortedByCoord.out.bam//')
        echo "Processing $prefix $sample"

        featureCounts -T ${nThreads} -a ${gtfFile} -o ${targetDir}/${sample}_raw_counts.txt $f
    done
done

