## Start with the genotype calls from Jacqueline Robinson

## Get separate set of calls for Amboseli & unadmixed individuals
sbatch --array=1-20 --mem=12G run.01.get_unadmixed.sh
sbatch --array=1-20 --mem=12G run.01.get_amboseli.sh

## Get set of mutual calls, and merged file
sbatch --array=1-20 --mem=12G run.02.merge.sh

## merge vcfs for yellow and anubis baboons
module load bcftools
bcftools concat /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed/02.anu.*.vcf.gz -O z -o /data/tunglab/tpv/panubis1_genotypes/anubis.vcf.gz
bcftools concat /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed/02.yel.*.vcf.gz -O z -o /data/tunglab/tpv/panubis1_genotypes/yellow.vcf.gz

## Estimate relatedness for unadmixed individuals using KING (http://people.virginia.edu/~wc9c/KING/manual.html)
# convert to plink format
plink --vcf hicov.amboseli.recode.vcf --maf 0.05 --recode --out plink.ambo 

# Install KING
##### wget http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz; tar -xvf Linux-king.tar.gz
# run KING
/data/tunglab/tpv/Programs/KING/king
module load gcc; module load glibc/2.18-gcb01

## Let's get just a bit of individual and site level information

module load vcftools
vcftools --gzvcf ./anubis.vcf.gz --depth --out anubis
vcftools --gzvcf ./yellow.vcf.gz --depth --out yellow
