## Get recombination maps for the panubis1 genome
# We'll be making 4 maps:
#        -anubis map (n=24 SW anubis)
#        -yellow map (n=XX high-cov yellow)
#        -SW yellow map
#        -Mikumi yellow map
#        -Downsampled anubis maps (to match partial yellow maps)
# Uses genotype calls for anubis and yellow individuals from /Amboseli_LocalAncestry/Resouces/Genotype_Calls

cd /data/tunglab/tpv/panubis1_genotypes/recombination

## Starting with the high coverage samples (in /data/tunglab/tpv/panubis1_genotypes/high_coverage and on github Amboseli_LocalAncestry/Resources/High Coverage Samples/)
module load bcftools; bcftools concat /data/tunglab/tpv/panubis1_genotypes/high_coverage/02.unadmixed_highcov.*.vcf.gz -O z -o ./unadmixed_highcov.vcf.gz

## Recombination map-specific genotype filtering for each dataset
# Filtering instructions from Pfeifer 2020 (vervet monkey genetic map)

## anubisSW
# Excluded sites with missing data to avoid errors and biases resulting from computational imputation
module load vcftools; vcftools --gzvcf unadmixed_highcov.vcf.gz --keep ../00_anu_SW.list --max-missing 1 --max-alleles 2 --recode --recode-INFO-all --out anubisSW --maf 0.05
## Remove clusters of SNPs within a 10bp window (GATK clusterSize=3, clusterWindowSize=10) and singletons within the target population
module load tabix; bgzip anubisSW.recode.vcf; module load java/1.8.0_45-fasrc01; module load tabix; tabix ./anubisSW.recode.vcf.gz
java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ../Panubis1.nochromname.fa -V ./anubisSW.recode.vcf.gz -cluster 3 -window 10 -o ./anubisSW_cluster.vcf.gz
# Stringent filtering for false variants
## SNPs showing an excess of heterozygosity were removed. Specifically, a P-value for Hardy–Weinberg Equilibrium was calculated using the “—hardy” option in VCFtools v.0.1.13 (Danecek et al. 2011), and SNPs with P < 0.01 removed.
module load vcftools; vcftools --hardy --gzvcf anubisSW_cluster.vcf.gz --out anubisSW

## SNPs that could be recipricolly lifted over with the human genome? We'll ignore this, assuming the quality of the baboon genome is sufficient. 

# Remove fixed alleles




## Phasing with beagle
############ Should we ask about error/differences in imputation based on the sample size (i.e. full panel vs partial yellow panels?) I'm actually going to just assume that's accounted for in the downsampling... 



