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

## Let's get just a bit of individual and site level information

module load vcftools
vcftools --gzvcf ./anubis.vcf.gz --depth --out anubis
vcftools --gzvcf ./yellow.vcf.gz --depth --out yellow

## for reference allele frequencies, see https://github.com/TaurVil/Amboseli_LocalAncestry/tree/master/Resources/unadmixed_reference_panels
