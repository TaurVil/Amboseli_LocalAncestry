## Goal: create vcf files of anubis and yellow baboons after masking 

## Step 1: get tracts to mask per individual 
mkdir mask_tracts_by_sample; mkdir indiv_vcfs 
for samp in `cat 00_refpanel.names`; do for version in `cat 00_versions.list`; do grep $samp ./masking_tracts/$version.bed > mask_tracts_by_sample/to_remove.$samp.$version.bed; done; echo $samp; done 
## Replace empty files with a minimal file (just excluding snps in the first BP). Doesn't delete information, but makes the files accessible to bedtools/vcftools down the road because they're still in the right format
for f in `find . -name '*.bed' -size 0`; do echo -e "chr1\t1\t2" > $f; done 

## Step 2: get individual vcf files 
## this is a lot of jobs, so to run it I actually split `cat 00_versions.list` into `head -N 00_versions.list | tail -5` to run sets of 5 (N=5,10,15,20)
for samp in `cat 00_refpanel.names`; do for mask in `head -20 00_versions.list | tail -5`; do sed -e s/SAMPLE_NAME/$samp/g -e s/MASK_NAME/$mask/g ./run.01.get_indiv_vcf.sh > ./r.$samp.$mask.sh; sbatch --mem=3G --nice ./r.$samp.$mask.sh; rm ./r.$samp.$mask.sh; done; done   

## Step 3: merge individual masked vcfs into refpanel vcfs 
mkdir refpanel_vcfs; for vers in `cat 00_versions.list`; do sed -e s/VERSION_NAME/$vers/g run.02.merge_filter_vcfs.sh > r.$vers.sh; sbatch --mem=8G r.$vers.sh; rm r.$vers.sh; done 


## Cleanup:
# We'll keep the masked refpanel vcfs, but remove the individual vcfs and the individual tracts to remove. 
rm -r mask_tracts_by_sample/
rm -r indiv_vcfs/
# Remove pre-filtered vcf files 
rm refpanel_vcfs/masked.yes*vcf.gz ; rm refpanel_vcfs/masked.no*vcf.gz
# also remove the record of singleton variants and variants less with less than 10 samples in one of the species: 
rm 01.*.to_remove.txt
