#!/bin/bash
cores=24
mem=100

new_sample_dir=
transcriptome_rna_cr_ref=
atac_cr_ref=

# non treated scRNA seq
sname="105_2g_DMSO"
cellranger count --id=${sname} \
    --transcriptome=${transcriptome_rna_cr_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem \
    --include-introns=true \
    --r1-length 26

sname="105_2g_DMSO_sub_0_4"
cellranger count --id=${sname} \
    --transcriptome=${transcriptome_rna_cr_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem \
    --include-introns=true \
    --r1-length 26


sname="HD_7HA_DMSO"
cellranger count --id=${sname} \
    --transcriptome=${transcriptome_rna_cr_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem \
    --include-introns=true \
    --r1-length 26

# senolyitic treated scRNA seq
sname="105_2g_ABT_263"
cellranger count --id=${sname} \
    --transcriptome=${transcriptome_rna_cr_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem \
    --include-introns=true \
    --r1-length 26

sname="HD_7HA_ABT_263"
cellranger count --id=${sname} \
    --transcriptome=${transcriptome_rna_cr_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem \
    --include-introns=true \
    --r1-length 26

# non treated snATACseq
sname="105_2g"
cellranger-atac count --id=${sname} \
    --reference=${atac_cr_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem 

sname="HD_7HA"
cellranger-atac count --id=${sname} \
    --reference=${atac_cr_ref} \
    --fastqs=${new_sample_dir} \
    --sample=${sname} \
    --localcores=$cores \
    --localmem=$mem 
