#!/bin/bash


source_folder=
genomePath="STAR_2710"
outputDir="${source_folder}/05.Star_mapped"
sourceDir="${source_folder}/04.Subsampled"

if [ ! -d $outputDir ]; then
    mkdir -p $outputDir
fi

for sample in $(ls $sourceDir); do
    if [ ! -d $outputDir/$sample ]; then
        mkdir -p $outputDir/$sample
    fi

    if [ ! -d $sourceDir/$sample ]; then
        continue
    fi

    if [ $sample == "FASTQC_sub_multiqc_data" ]; then
        continue
    fi

    for f in $(ls $sourceDir/$sample/*R1*gz); do
        prefix=$(echo $f | rev | cut -d / -f 1 | rev | sed 's/_R1_001.fastq.gz//')
        r1=$sourceDir/$sample/${prefix}_R1_001.fastq.gz
        r2=$sourceDir/$sample/${prefix}_R2_001.fastq.gz

        echo "Processing $prefix $sample"

        STAR 	--runThreadN 16 \
                --genomeDir $genomePath \
                --readFilesIn $r1 $r2 \
                --outSAMtype BAM SortedByCoordinate \
                --runMode alignReads \
                --outFileNamePrefix ${outputDir}/${sample}/${prefix}_ \
                --outReadsUnmapped Fastx \
                --readFilesCommand zcat   
    done
done
