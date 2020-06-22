#!/bin/bash

module load java; module load samtools

g=${SLURM_ARRAY_TASK_ID}

f=`head -$g 00_bams.txt | tail -1`;

path_genome=/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa

#samtools index ${f}_chr.bam

samtools view -b -L 01_targetted_chroms.bed ${f}_chr.bam > target.$f.bam
samtools index target.$f.bam

module load GATK

GenomeAnalysisTK.sh HaplotypeCaller -ERC GVCF -I target.$f.bam  -R $path_genome  -O gVCF/$f.g.vcf.gz -mbq 20

