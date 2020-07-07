## Goal: get genotypes for high coverage, unadmixed individuals (both yellow, anubis, and other)
# Used for recombiantion maps and ibdmix

mkdir /data/tunglab/tpv/panubis1_genotypes/high_coverage; cd /data/tunglab/tpv/panubis1_genotypes/high_coverage

sbatch --array=1-20 --mem=6G run.get_vcf.sh


## Get all non-singleton sites covered in all high coverage samples and filtered based on GATK calls and to remove chunks of nearby variants
