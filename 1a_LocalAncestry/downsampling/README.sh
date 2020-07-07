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
## merge with refpanel and get genolik files
sbatch --array=1-6 --mem=24G run.05.merge_get_genolik.sh

## use vcftools to get the depth (idepth), relatedness, and individual names from one of hte CommonCalls vcfs. These aren't necessary for anything, but help figure out which individuals are what for the *.h files

## 32 to 38 are the amboseli samples
## SW are up to 31, mixed between yellow (n=7) and anubis (n=24)
## 7 Mara: 39 to 45
## 15 Mikumi: 46 to 59
## 6 Tulane: 60 to 65
## 4 BGP anubis: 66 to 69
## 1 BGP yellow: 70 (71 is a copy of Mikumi 7)

## Get LCLAE calls
mkdir raw_calls
for cov in `cat 01_sets.txt`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 71 02.anubis.h 02.yellow.h $id | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/full_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; echo $id $name;  done; echo $cov; done 

## add chromosome names, and merge into one file
for cov in `cat 01_sets.txt`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do sed -i 's/^/\t/' ./raw_calls/full_refpanel.35kb.p2.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/full_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; cat ./raw_calls/full_refpanel.35kb.p2.$name.$cov.chr* > ./raw_calls/full_refpanel.35kb.p2.$name.$cov.txt; rm ./raw_calls/full_refpanel.35kb.p2.$name.$cov.chr*txt; done; done

## Get majority rule and tracts
module load R; for name in `cat 00_names.txt`; do for covs in `cat 01_sets.txt`; do sed -e s/INDIV/$name/g -e s/COVERAGE/$covs/g run.07.calls_to_MajRule_tracts.R > r.$covs.R; sbatch r.$covs.R; rm r.$covs.R; done; done


## get calls for sw and wall refpanels
for cov in `echo original 2x 1x`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 71 02.anubis_sw.h 02.yellow_sw.h $id | ./geno_lik2 .2 35000 > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; echo $id $name;  done; echo $cov; done 



for cov in `echo original 2x 1x`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 71 02.anubis_sw.h 02.yellow_sw.h $id | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; echo $name $cov $chrom ; sed -i 's/^/\t/' ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; echo $id $name; cat ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.chr* > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.txt; rm ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.chr*txt; done; echo $cov; done 

c1-01-4, but not c1-01-3? 


for cov in `echo 2x 1x`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 71 02.anubis_sw.h 02.yellow_sw.h $id | ./lclae_master/source/geno_lik2 .2 35000 > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; echo $chrom; done; echo $id $name;  done; echo $cov; done

for cov in `echo 2x 1x`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 71 02.anubis_wall.h 02.yellow_wall.h $id | ./lclae_master/source/geno_lik2 .2 35000 > ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; echo $chrom; done; echo $id $name;  done; echo $cov; done

## add chrom and merge
for cov in `echo 2x 1x original`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do sed -i 's/^/\t/' ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; cat ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.chr* > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.txt; rm ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.chr*txt; done; done

for cov in `echo 2x 1x original`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do sed -i 's/^/\t/' ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; cat ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.chr* > ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.txt; rm ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.chr*txt; done; done

# correct ones that were run twice and therefore duplicated: for f in `ls ./raw_calls/wall_refpanel*`; do sed -i s/chr.*chr/chr/g $f; done

# get majority rule and tracts
sed -e s/full_refpanel/sw_refpanel/g run.07.calls_to_MajRule_tracts.R > run.07a.sw_refpanel.R
sed -e s/full_refpanel/wall_refpanel/g run.07.calls_to_MajRule_tracts.R > run.07b.wall_refpanel.R
module load R
module load R; for name in `cat 00_names.txt`; do for covs in `cat 01_sets.txt`; do sed -e s/INDIV/$name/g -e s/COVERAGE/$covs/g run.07a.sw_refpanel.R > r.$covs.R; sbatch r.$covs.R; rm r.$covs.R; done; done
module load R; for name in `cat 00_names.txt`; do for covs in `cat 01_sets.txt`; do sed -e s/INDIV/$name/g -e s/COVERAGE/$covs/g run.07b.wall_refpanel.R > r.$covs.R; sbatch r.$covs.R; rm r.$covs.R; done; done
