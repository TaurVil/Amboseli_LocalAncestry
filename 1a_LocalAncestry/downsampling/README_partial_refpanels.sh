## get ancestry calls for sw and wall refpanels
for cov in `echo original 2x 1x`; do 
  for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; 
    for chrom in `cut -f 1 01_targetted_chroms.bed`; do 
    touch ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; l=`wc -l ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt | cut -f 2 | sed 's/ .*//g'`; if (($l==0)); then echo "STARTING"; 
    grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 70 02.anubis_sw.h 02.yellow_sw.h $id | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; fi; echo "done" $chrom; done; echo $id $name;  done; echo $cov; done 

for cov in `echo original 2x 1x`; do 
  for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; 
    for chrom in `cut -f 1 01_targetted_chroms.bed`; do 
    touch ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; l=`wc -l ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt | cut -f 2 | sed 's/ .*//g'`; if (($l==0)); then echo "STARTING"; 
    grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 70 02.anubis_wall.h 02.yellow_wall.h $id | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; fi; echo "done" $chrom; done; echo $id $name;  done; echo $cov; done

for cov in `echo original 2x 1x`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 71 02.anubis_sw.h 02.yellow_sw.h $id | ./lclae_master/source/geno_lik2 .2 35000 > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; echo $chrom; done; echo $id $name;  done; echo $cov; done
for cov in `echo original 2x 1x`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do grep ^$chrom genolik.$cov.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 71 02.anubis_wall.h 02.yellow_wall.h $id | ./lclae_master/source/geno_lik2 .2 35000 > ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; echo $chrom; done; echo $id $name;  done; echo $cov; done

## add chrom and merge
for cov in `echo 2x 1x original`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do sed -i 's/^/\t/' ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; cat ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.chr* > ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.txt; rm ./raw_calls/sw_refpanel.35kb.p2.$name.$cov.chr*txt; done; done
for cov in `echo 2x 1x original`; do for id in `seq 32 38`; do name=`head -$id tmp.012.indv | tail -1 | sed -e s/.varia.*$//`; for chrom in `cut -f 1 01_targetted_chroms.bed`; do sed -i 's/^/\t/' ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.$chrom.txt; done; cat ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.chr* > ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.txt; rm ./raw_calls/wall_refpanel.35kb.p2.$name.$cov.chr*txt; done; done

# correct ones that were run twice and therefore duplicated: for f in `ls ./raw_calls/wall_refpanel*`; do sed -i s/chr.*chr/chr/g $f; done

# get majority rule and tracts
sed -e s/full_refpanel/sw_refpanel/g run.07.calls_to_MajRule_tracts.R > run.07a.sw_refpanel.R
sed -e s/full_refpanel/wall_refpanel/g run.07.calls_to_MajRule_tracts.R > run.07b.wall_refpanel.R
module load R
module load R; for name in `cat 00_names.txt`; do for covs in `cat 01_sets.txt`; do sed -e s/INDIV/$name/g -e s/COVERAGE/$covs/g run.07a.sw_refpanel.R > r.$covs.R; sbatch --mem=8G r.$covs.R; rm r.$covs.R; done; done
module load R; for name in `cat 00_names.txt`; do for covs in `cat 01_sets.txt`; do sed -e s/INDIV/$name/g -e s/COVERAGE/$covs/g run.07b.wall_refpanel.R > r.$covs.R; sbatch --mem=8G r.$covs.R; rm r.$covs.R; done; done
