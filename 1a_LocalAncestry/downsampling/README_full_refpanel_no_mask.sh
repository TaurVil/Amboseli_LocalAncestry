## Working with 7 high coverage amboseli individuals, the bams for which are in /data/tunglab/tpv/mapped_bams/amboseli_hicov/

cd /data/tunglab/tpv/local_ancestry/downsampling

## Get coverage for each sample
## first command in this script. Comment the first one out. 
ls /data/tunglab/tpv/mapped_bams/amboseli_hicov/ | grep -v 'bai' | sed s/_MarkDuplicates.bam//g > 00_names.txt 
for f in `cat 00_names.txt`; do sed -e s/FILENAME/$f/g run.01.calculate_coverage.sh > tmp.sh; sbatch --mem=8G --out cov.$f.out tmp.sh; rm tmp.sh; done

## Downsample to get 10x data, as well as 5x, 2x, 1x, and p5x.
cat *out | sed 's/Average =  //g' > 00_covs.txt
sbatch --array=1-7 --mem=8G run.02.downsample.bam
## check the output for coverage information
## second command in this script. Comment the first one out. 
for f in `cat 00_bams.txt`; do sed -e s/FILENAME/$f/g run.01.calculate_coverage.sh > tmp.sh; sbatch --mem=8G --out cov.$f.out tmp.sh; rm tmp.sh; done

## Convert bams to chrXX rather than just XX
for file in AMB*.bam; do filename=`echo $file | sed s/.bam//g`; samtools view -H $file | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > ${filename}_chr.bam; echo $file; echo ${filename}_chr.bam; done

## Call genotypes for each set of bams
mkdir gVCF; sbatch --array=1-35 --mem=16G  run.03.gvcf.sh
sbatch --array=1-7 --mem=16G  run.03.gvcf_original.sh

## Joint genotype calling for each set of bams 
sbatch --array=1-6 --mem=16G run.04.merge_gvcfs.sh
## Clean out all pre-requisite files up to this point. 

# For the refpanel, get the list of kept sites, and adjust to use the 'chr' nomenclature
## in /data/tunglab/tpv/local_ancestry/unadmixed_individuals/refpanel.kept.sites
https://github.com/TaurVil/Amboseli_LocalAncestry/blob/master/Resources/Genotype_Calls/make_merged_refpanel.sh 

## merge with refpanel and get genolik files
sbatch --array=1-6 --mem=24G run.05.merge_get_genolik.sh

## use vcftools to get the depth (idepth), relatedness, and individual names from one of the CommonCalls vcfs. These aren't necessary for anything, but help figure out which individuals are what for the *.h files

## 32 to 38 are the amboseli samples
## SW are up to 31, mixed between yellow (n=7) and anubis (n=24)
## 7 Mara: 39 to 45
## 15 Mikumi: 46 to 59
## 6 Tulane: 60 to 65
## 4 BGP anubis: 66 to 69
## 1 BGP yellow: 70

## Get LCLAE calls
mkdir raw_calls
for cov in `cat 01_sets.txt`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do echo "starting" $chrom; grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 70 02.anubis.h 02.yellow.h $id | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/full_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; echo $id $name;  done; echo $cov; done 

## add chromosome names, and merge into one file
for cov in `cat 01_sets.txt`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do sed -i 's/^/\t/' ./raw_calls/full_refpanel.35kb.p2.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/full_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; cat ./raw_calls/full_refpanel.35kb.p2.$name.$cov.chr* > ./raw_calls/full_refpanel.35kb.p2.$name.$cov.txt; rm ./raw_calls/full_refpanel.35kb.p2.$name.$cov.chr*txt; done; done

## Get majority rule and tracts
module load R; for name in `cat 00_names.txt`; do for covs in `cat 01_sets.txt`; do sed -e s/INDIV/$name/g -e s/COVERAGE/$covs/g run.07.calls_to_MajRule_tracts.R > r.$covs.R; sbatch --mem=8G r.$covs.R; rm r.$covs.R; done; done
