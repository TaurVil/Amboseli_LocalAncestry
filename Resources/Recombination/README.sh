## Get recombination maps for the panubis1 genome
# We'll be making 4 maps:
#        -anubis map (n=24 SW anubis)
#        -yellow map (n=XX high-cov yellow)
#        -SW yellow map
#        -Mikumi yellow map
#        -Downsampled anubis maps (to match partial yellow maps)
# Uses genotype calls for anubis and yellow individuals from /Amboseli_LocalAncestry/Resouces/Genotype_Calls

pwd /data/tunglab/tpv/panubis1_genotypes/recombination

## Filter genotypes
# Filtering instructions from Pfeifer 2020 (vervet monkey genetic map)
## Remove clusters of SNPs within a 10bp window (GATK clusterSize=3, clusterWindowSize=10) 

module load java/1.8.0_45-fasrc01; module load tabix
tabix ../anubis.vcf.gz
java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ../Panubis1.nochromname.fa -V ../anubis.vcf.gz --cluster-size 3 --cluster-window-size 10 -o ./anubis_cluster.vcf.gz



# Excluded sites with missing data (n=XX, XX%) to avoid errors and biases resulting from computational imputation
# Stringent filtering for false variants
## SNPs showing an excess of heterozygosity were removed. Specifically, a P-value for Hardy–Weinberg Equilibrium was calculated using the “—hardy” option in VCFtools v.0.1.13 (Danecek et al. 2011), and SNPs with P < 0.01 removed.

## SNPs that could be recipricolly lifted over with the human genome?

# Remove fixed alleles




## Phasing with beagle
############ Should we ask about error/differences in imputation based on the sample size (i.e. full panel vs partial yellow panels?) I'm actually going to just assume that's accounted for in the downsampling... 



