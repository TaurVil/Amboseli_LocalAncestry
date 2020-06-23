#!/bin/bash

index=${SLURM_ARRAY_TASK_ID}

coverage=`head -$index 01_sets.txt | tail -1`
path_genome=/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa

 module load vcftools; module load tabix

vcftools --gzvcf ./test.$coverage.vcf.gz --max-alleles 2 --minQ 30 --positions ../refpanel.kept.sites --recode --out filt.$coverage
bgzip filt.$coverage.recode.vcf
tabix ./filt.$coverage.recode.vcf.gz

module load java/1.8.0_45-fasrc01

java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T CombineVariants -R $path_genome --variant ./filt.$coverage.recode.vcf.gz --variant ../refpanel.vcf.gz -o merged.$coverage.vcf.gz -genotypeMergeOptions UNIQUIFY -L 01_targetted_chroms.bed

java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T SelectVariants -R $path_genome -V merged.$coverage.vcf.gz -select 'set == "Intersection"' -o CommonCalls.$coverage.vcf

vcftools --vcf CommonCalls.$coverage.vcf --max-alleles 2 --recode --out CommonCalls.biallelic.$coverage

module load bcftools
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL CommonCalls.biallelic.$coverage.recode.vcf | sed 's/|/\//g' > merged.$coverage.forgenolik.vcf

sed '/^#/d' merged.$coverage.forgenolik.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /data/tunglab/tpv/Programs/LCLAE/filtbaboon1b 71 > genolik.$coverage.genolik

