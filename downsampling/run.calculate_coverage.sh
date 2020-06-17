#!/bin/bash

module load samtools

name=FILENAME

samtools depth /data/tunglab/tpv/mapped_bams/amboseli_hicov/$name*.bam  > cov.$name.bam ; awk '{sum+=$3} END { print "Average = ",sum/NR}' cov.$name.bam; rm cov.$name.bam
