#!/bin/bash
cores=24
mem=100

new_sample_dir=
cellranger_atac_ref=

sname="105_2g_1"
echo $sname
/home/am3019/software/cellranger-atac-2.1.0/cellranger-atac count --id=${sname} \
    --reference=${cellranger_atac_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem 

sname="105_2g_2"
echo $sname
/home/am3019/software/cellranger-atac-2.1.0/cellranger-atac count --id=${sname} \
    --reference=${cellranger_atac_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem 



sname="HD_74A_1"
echo $sname
/home/am3019/software/cellranger-atac-2.1.0/cellranger-atac count --id=${sname} \
    --reference=${cellranger_atac_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem 

sname="HD_74A_2"
echo $sname
/home/am3019/software/cellranger-atac-2.1.0/cellranger-atac count --id=${sname} \
    --reference=${cellranger_atac_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem 



