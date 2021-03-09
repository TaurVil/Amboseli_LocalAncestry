## Run Ancestry HMM for all Amboseli individuals

cd /data/tunglab/tpv/local_ancestry/HMM

### Prep data from vcf file 
cd /data/tunglab/tpv/local_ancestry/HMM/
module load bcftools; bcftools query -f '%CHROM \t %POS \t %REF \t %ALT [\t %AD]\n' /data/tunglab/tpv/panubis1_genotypes/masked_final/final.amboseli_with_masked_refpanel.vcf.gz > temp
sed -i 's/\./0,0/g' temp
sed -i 's/,/\t/g' temp

module load bcftools ; bcftools query -l /data/tunglab/tpv/panubis1_genotypes/masked_final/final.amboseli_with_masked_refpanel.vcf.gz > samples.list


## Prepare to run Ancestry HMM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/tunglab/tpv/Programs/conda/lib/
export PATH=$PATH:/data/tunglab/tpv/Programs/conda/lib/
module load Anaconda; conda --version
module load libtool; libtool --version

# Install the AncestryHMM software
# wget https://github.com/russcd/Ancestry_HMM/archive/master.zip; gunzip master.zip
# cd Ancestry_HMM/src/
# make






# Create ancestry tracts from posterior calls
get_ancestryHMM_tracts.R with INDIV and CHROM 


# ran using array_for_get_ancestry.sh, uses 00_ambo_names.list and the chromosome-specific version of get_ancestryHMM_tracts.R (get.CHROM.R)
sbatch --array=1-445 array_for_get_ancestry.sh
# outputs a file for CHROM_INDIV in tracts

# Convert tracts to HMM tracts
for i in `cat 00_ambo_names.list`; do cat tracts/*${i}* > HMMtracts/${i}.txt; done 

## We'll pull these to the local directory to analyze
