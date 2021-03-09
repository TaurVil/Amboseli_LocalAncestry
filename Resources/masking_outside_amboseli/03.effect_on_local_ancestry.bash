## Goal: test the impacts of masking on calling local ancestry

## we'll just run this on chroms 17-20, so we'll be using the `targets` version of the refpanel vcfs

## Step 1: get genolik files for LCLAE
## merge with genotype information for downsampled Amboseli samples (1x, 2x, 5x, and original coverage) 
## Creates genolik files. We'll remove everything except those after this runs successfully. 
mkdir merged_to_call
for vers in `cat 00_versions.list`; do sed -e s/VERSION_NAME/$vers/g run.03.merge_get_genolik.sh > r.$vers.sh; sbatch --array=1-4 --mem=16G r.$vers.sh; rm r.$vers.sh; done 

## Step 2: get chromosome-level genolik files 
## Get genolik by chromosome, why are some of these empty? (turned out to be .: genotypes rather than ./.:)
for chrom in `cat ../downsampling/01_targetted_chroms.bed | cut -f 1`; do for vers in `cat 00_versions.list`; do for cov in `cat 00_coverages.list`; do grep chr${chrom} merged_to_call/genolik.$cov.$vers.genolik > genolik.$vers.$cov.$chrom.genolik; done; done; done 

## Step 3: get raw calls
mkdir raw_calls 
for vers in `cat 00_versions.list`; do for cov in `cat 00_coverages.list`; do sed -e s/VERSION_NAME/$vers/g -e s/COVERAGE_NAME/$cov/g run.04.get_raw_calls.sh > r.$cov.$vers.sh; sbatch --array=32-38 --nice r.$cov.$vers.sh; rm r.$cov.$vers.sh; done; done 
## iterate until all run successfully 

## Step 4: convert raw calls into tracts 
### add chromosome and merge for each individual/condition/coverage 
for cov in `cat 00_coverages.list`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for version in `cat 00_versions.list`; do for chrom in `seq 17 20`; do sed -i 's/^/\t/' ./raw_calls/full_refpanel.35kb.p2.$id.$chrom.$cov.$version.txt; sed -i s/^/chr${chrom}/ ./raw_calls/full_refpanel.35kb.p2.$id.$chrom.$cov.$version.txt; done; cat ./raw_calls/full_refpanel.35kb.p2.$id.*.$cov.$version.txt > ./raw_calls/full_refpanel.35kb.p2.$name.$cov.$version.txt; mv ./raw_calls/full_refpanel.35kb.p2.$id.*.$cov.$version.txt used/; echo $name $id $cov $version; done; done; done 
mkdir tracts 
module load R; for covs in `cat 00_coverages.list`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for version in `cat 00_versions.list`; do sed -e s/INDIV/$name/g -e s/VERSION/$version/g -e s/COVERAGE/$covs/g run.07.calls_to_MajRule_tracts.R > r.$covs.R; sbatch --mem=8G  --nice r.$covs.R; rm r.$covs.R; done; done; done 


## Cleanup
# remove raw LCLAE calls
rm -r raw_calls/
# zip up the genolik files 
tar -czvf genoliks.tar.gz merged_to_call/ ; rm -r merged_to_call/

