#!/bin/bash

module load java; module load samtools

g=${SLURM_ARRAY_TASK_ID}

f=`head -$g 00_names.txt | tail -1`;

path_genome=/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa

file=/data/tunglab/tpv/mapped_bams/amboseli_hicov/$f*.bam

#samtools view -H $file | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > $f.original.bam

#samtools index $f.original.bam


samtools view -b -L 01_targetted_chroms.bed $f.original.bam > target.$f.bam
samtools index target.$f.bam


module load GATK

GenomeAnalysisTK.sh HaplotypeCaller -ERC GVCF -I target.$f.bam -R $path_genome  -O gVCF/$f.original.g.vcf.gz -mbq 20

