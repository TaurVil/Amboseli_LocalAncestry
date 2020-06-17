#!/bin/bash

g=${SLURM_ARRAY_TASK_ID}

name=`head -$g 00_names.txt | tail -1`
i_cov=`head -$g 00_covs.txt | tail -1`


perc=`echo 10/$i_cov | bc -l`

module load samtools
samtools view -s $perc -b /data/tunglab/tpv/mapped_bams/amboseli_hicov/$name*.bam > $name.10x.bam

perc2=`echo $perc/2 | bc -l`
samtools view -s $perc2 -b /data/tunglab/tpv/mapped_bams/amboseli_hicov/$name*.bam > $name.5x.bam


perc2=`echo $perc/5 | bc -l`
samtools view -s $perc2 -b /data/tunglab/tpv/mapped_bams/amboseli_hicov/$name*.bam > $name.2x.bam

perc2=`echo $perc/10 | bc -l`
samtools view -s $perc2 -b /data/tunglab/tpv/mapped_bams/amboseli_hicov/$name*.bam > $name.1x.bam

perc2=`echo $perc/20 | bc -l`
samtools view -s $perc2 -b /data/tunglab/tpv/mapped_bams/amboseli_hicov/$name*.bam > $name.p5x.bam

