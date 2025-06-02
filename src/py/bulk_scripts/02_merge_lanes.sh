#!/bin/bash

projectFolder=
sourceFolder=${projectFolder}/02.Renamed
outputFolder=${projectFolder}/03.LanesMerged

for sampleName in $(ls ${sourceFolder}); do
    if [ ! -d ${sourceFolder}/${sampleName} ]; then
        continue
    fi

    if [ ! -d ${outputFolder}/${sampleName} ]; then
        mkdir -p ${outputFolder}/${sampleName}
    fi

    if [ ${sampleName} == "FASTQC_multiqc_data" ]; then
        continue
    fi

    echo "Processing ${sampleName}"
    cp ${sourceFolder}/${sampleName}/*.fastq.gz ${outputFolder}/${sampleName}

    if [ -f ${sourceFolder}/${sampleName}/${sampleName}_L002_R1_001.fastq.gz ]; then
        pigz -d -p 24 ${outputFolder}/${sampleName}/${sampleName}_L001_R1_001.fastq.gz
        pigz -d -p 24 ${outputFolder}/${sampleName}/${sampleName}_L001_R2_001.fastq.gz
        pigz -d -p 24 ${outputFolder}/${sampleName}/${sampleName}_L002_R1_001.fastq.gz
        pigz -d -p 24 ${outputFolder}/${sampleName}/${sampleName}_L002_R2_001.fastq.gz

        echo $(wc -l ${outputFolder}/${sampleName}/${sampleName}_L001_R1_001.fastq)
        
        cat ${outputFolder}/${sampleName}/${sampleName}_L001_R1_001.fastq ${outputFolder}/${sampleName}/${sampleName}_L002_R1_001.fastq > temp.fastq
        mv temp.fastq ${outputFolder}/${sampleName}/${sampleName}_L001_R1_001.fastq
        rm ${outputFolder}/${sampleName}/${sampleName}_L002_R1_001.fastq

        cat ${outputFolder}/${sampleName}/${sampleName}_L001_R2_001.fastq ${outputFolder}/${sampleName}/${sampleName}_L002_R2_001.fastq > temp.fastq
        mv temp.fastq ${outputFolder}/${sampleName}/${sampleName}_L001_R2_001.fastq
        rm ${outputFolder}/${sampleName}/${sampleName}_L002_R2_001.fastq

        echo $(wc -l ${outputFolder}/${sampleName}/${sampleName}_L001_R1_001.fastq)

        pigz -p 24 ${outputFolder}/${sampleName}/${sampleName}_L001_R1_001.fastq
        pigz -p 24 ${outputFolder}/${sampleName}/${sampleName}_L001_R2_001.fastq
    fi

done