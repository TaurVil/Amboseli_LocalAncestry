#!/bin/bash

index=${SLURM_ARRAY_TASK_ID}
coverage=`head -$index 01_sets.txt | tail -1`
path_genome=/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa

module load java; module load samtools; module load python
module load htslib
module load GATK

ls gVCF/*.$coverage*vcf.gz > 02.$coverage.list


for chrom in `cut -f 1 01_targetted_chroms.bed`; do GenomeAnalysisTK.sh CombineGVCFs  -R $path_genome -O test.$coverage.$chrom.g.vcf.gz -V 02.$coverage.list -L $chrom; GenomeAnalysisTK.sh GenotypeGVCFs -R $path_genome -V test.$coverage.$chrom.g.vcf.gz -O test.$coverage.$chrom.vcf.gz --tmp-dir /data/tunglab/tpv/scratch/; rm test.$coverage.$chrom.g.vcf.gz; rm test.$coverage.$chrom.g.vcf.gz.tbi; done

rm 02.$coverage.list; 
