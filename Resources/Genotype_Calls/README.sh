## Start with the genotype calls from Jacqueline Robinson

## Get separate set of calls for Amboseli & unadmixed individuals
sbatch --array=1-20 --mem=12G run.01.get_unadmixed.sh
sbatch --array=1-20 --mem=12G run.01.get_amboseli.sh

## Get set of mutual calls, and merged file
sbatch --array=1-20 --mem=12G run.02.merge.sh


## Let's get just a bit of individual and site level information

module load vcftools
vcftools --gzvcf calls_unadmixed/01b.gatk_filtered.unadmixed.20.vcf.gz --depth --out chr20_unadmixed
