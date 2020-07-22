#!/bin/bash
path_genome=/data/tunglab/shared/genomes/panubis1/Panubis_1.0.fa

## This will work as the chromosomes are indexed 1-20,X,Y,scafoldZZZ rather than "chr1", "chr2", etc.
index=${SLURM_ARRAY_TASK_ID}
chrom=$index

module load java/1.8.0_45-fasrc01; module load tabix
module load samtools; module load vcftools

# Get genotype calls for the high coverage samples, putatively unadmixed, samples.
vcftools --gzvcf /data/tunglab/asf40/wgs_data/from_Jacqueline/baboon1k_v1_snpEff_chr$chrom'.vcf.gz' --keep 10.hicov.txt --max-alleles 2 --remove-indels --recode --out ./01.$chrom.unadmixed --recode-INFO-all

# vcftools filtering pass 1
vcftools --gzvcf ./01.$chrom.unadmixed.recode.vcf --max-missing 1 --minQ 30 --max-alleles 2 --remove-indels --maf 0.001 --recode --out ./01b.$chrom.unadmixed --recode-INFO-all


gunzip ./01b.$chrom.unadmixed.recode.vcf
# GATK filtering for true variants and to remove clusters
java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ../Panubis1.nochromname.fa -V ./01b.$chrom.unadmixed.recode.vcf -filterName "FS" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -cluster 3 -window 10 -o ./01c.$chrom.unadmixed.vcf.gz

# Identify singletons
vcftools --gzvcf ./01c.$chrom.unadmixed.vcf.gz --singletons --out ./01c.$chrom.unadmixed

# Output filtered genotypes
vcftools --gzvcf ./01c.$chrom.unadmixed.vcf.gz --exclude-positions ./01c.$chrom.unadmixed.singletons --recode --out 02.unadmixed_highcov.$chrom --recode-INFO-all --remove-filtered-all

bgzip 02.unadmixed_highcov.$chrom.recode.vcf ; tabix 02.unadmixed_highcov.$chrom.recode.vcf.gz

bgzip ./01b.$chrom.unadmixed.recode.vcf; bgzip ./01.$chrom.unadmixed.recode.vcf
