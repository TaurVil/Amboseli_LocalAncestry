### This folder houses files of individuals from each possible reference population. 

## Anubis: n=41
# SW: n=24 from SW NPRC founder individuals Robinson et al. 2019
# Wall: n=13
# Wall, Masai Mara: n=7 from Wall et al. 2016
# Wall, Tulane NPRC: n=6 from Wall et al. 2016
# Baboon Genome Project: n=4 from Rogers et al. 2019

## Yellow: n=23
# SW: n=7 from SW NPRC founder individuals in Robinson et al. 2019. 
# Wall: n=15
# Mikumi: n=14, 10 high coverage new sequencing & 4 low coverage from Wall et al. 2016
# Baboon Genome Project: n=1 from Mikumi and Rogers et al. 2019
######## We are not using pcyn_16098, who is the same individual as Mik_007



## Get reference panel allele frequencies
cd /data/tunglab/tpv/panubis1_genotypes/
module load vcftools
# all anubis/yellow samples
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.anu.$chrom.recode.vcf.gz --keep 00_anu.list --freq --out ./refpanel_allele_freqs/all_anubis.$chrom ; echo $chrom; done
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.yel.$chrom.recode.vcf.gz --keep 00_yel.list --freq --out ./refpanel_allele_freqs/all_yellow.$chrom ; echo $chrom; done

# subsets of anubis individuals
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.anu.$chrom.recode.vcf.gz --keep 00_anu_SW.list --freq --out ./refpanel_allele_freqs/sw_anubis.$chrom ; echo $chrom; done
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.anu.$chrom.recode.vcf.gz --keep 00_anu_mara.list --freq --out ./refpanel_allele_freqs/mara_anubis.$chrom ; echo $chrom; done
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.anu.$chrom.recode.vcf.gz --keep 00_anu_tulane.list --freq --out ./refpanel_allele_freqs/tulane_anubis.$chrom ; echo $chrom; done
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.anu.$chrom.recode.vcf.gz --keep 00_anu_wall.list --freq --out ./refpanel_allele_freqs/wall_anubis.$chrom ; echo $chrom; done

# subsets of yellow individuals
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.yel.$chrom.recode.vcf.gz --keep 00_yel_SW.list --freq --out ./refpanel_allele_freqs/sw_yellow.$chrom ; echo $chrom; done
for chrom in `seq 1 20`; do vcftools --gzvcf ./calls_unadmixed/02.yel.$chrom.recode.vcf.gz --keep 00_yel_wall.list --freq --out ./refpanel_allele_freqs/wall_yellow.$chrom ; echo $chrom; done

## Remove logs, combine files, split columns apart
rm refpanel_allele_freqs/*log
ls refpanel_allele_freqs/*.9.frq > 10.ref_sets.list; sed -i s/.9.frq//g 10.ref_sets.list
for f in `ls refpanel_allele_freqs/*`; do sed -i '1d' $f; done
for f in `ls refpanel_allele_freqs/*`; do sed -i 's/:/\t/g' $f; done
for f in `cat 10.ref_sets.list`; do cat $f.* > $f.freq; done 

module load R; R
library(data.table)
fread("./refpanel_allele_freqs/wall_anubis.freq") -> wall_anubis
fread("./refpanel_allele_freqs/wall_yellow.freq") -> wall_yellow

fread("./refpanel_allele_freqs/tulane_anubis.freq") -> tulane_anubis
fread("./refpanel_allele_freqs/mara_anubis.freq") -> mara_anubis
fread("./refpanel_allele_freqs/sw_anubis.freq") -> sw_anubis
fread("./refpanel_allele_freqs/sw_yellow.freq") -> sw_yellow

fread("./refpanel_allele_freqs/all_anubis.freq") -> all_anubis
fread("./refpanel_allele_freqs/all_yellow.freq") -> all_yellow

colnames(all_anubis) <- c("chrom", "snp", "num_alleles", "n_all_anubis", "ref", "all_anubis_ref", "alt", "all_anubis")
all_anubis[,-c(3,6)] -> all

colnames(all_yellow) <- c("chrom", "snp", "num_alleles", "n_all_yellow", "ref", "ref", "alt", "all_yellow")
colnames(wall_anubis) <- c("chrom", "snp", "num_alleles", "n_wall_anubis", "ref", "ref", "alt", "wall_anubis")
colnames(wall_yellow) <- c("chrom", "snp", "num_alleles", "n_wall_yellow", "ref", "ref", "alt", "wall_yellow")
colnames(tulane_anubis) <- c("chrom", "snp", "num_alleles", "n_tulane_anubis", "ref", "ref", "alt", "tulane_anubis")
colnames(mara_anubis) <- c("chrom", "snp", "num_alleles", "n_mara_anubis", "ref", "ref", "alt", "mara_anubis")
colnames(sw_anubis) <- c("chrom", "snp", "num_alleles", "n_sw_anubis", "ref", "ref", "alt", "sw_anubis")
colnames(sw_yellow) <- c("chrom", "snp", "num_alleles", "n_sw_yellow", "ref", "ref", "alt", "sw_yellow")

all[,c(1:2,4:5,6,3)] -> all
cbind(all, all_yellow[,c(8,4)], wall_anubis[,c(8,4)], wall_yellow[,c(8,4)], 
      tulane_anubis[,c(8,4)], mara_anubis[,c(8,4)], 
      sw_anubis[,c(8,4)], sw_yellow[,c(8,4)]) -> all
rm(all_yellow, all_anubis, wall_anubis, wall_yellow, tulane_anubis, mara_anubis, sw_anubis, sw_yellow)
save.image("./reference_allele_freq.RData")

