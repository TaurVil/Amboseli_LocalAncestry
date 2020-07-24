## Can we detect introgression in our reference panel samples? 

cd /data/tunglab/tpv/local_ancestry/unadmixed_individuals

# Get genotype likelihoods for unadmixed individuals
module load bcftools
for index in `seq 1 20`; do bcftools merge /data/tunglab/tpv/panubis1_genotypes/chrom_vcfs/02.anu.$index.recode.vcf.gz /data/tunglab/tpv/panubis1_genotypes/chrom_vcfs/02.yel.$index.recode.vcf.gz -O z -o $index.vcf.gz; echo $index; done

vcftools --gzvcf 20.vcf.gz --012 --out refpanel

sbatch --array=1-20 --mem=5G run.get_genolik.sh

## 41 anubis 1-41
### SW: 1-24, Mara 25-31, Tulane 32-37, BGP 38-41
## yellow 42-63
### SW: 42-48, Mikumi 49-63

mkdir raw_calls
echo 41 `seq 1 41` > 02.anubis.h; echo 22 `seq 42 63` > 02.yellow.h
echo 21 `seq 42 63` > 02.yellow_minus1.h; echo 40 `seq 1 41` > 02.anubis_minus1.h

sbatch --array=42-62 run.get_raw_calls.sh
sbatch --array=1-40 run.get_raw_calls2.sh

## 63 is special and the *.h file needs to be done manually
f=63
for chrom in `seq 1 20 `; do touch ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; l=`wc -l raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt | cut -f 2 | sed 's/ .*//g'`; if (($l==0)); then echo "STARTING"; cat genolik.$chrom.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 63 02.anubis.h $f.h $f | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; fi ; echo $chrom; done

f=41
for chrom in `seq 1 20 `; do touch ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; l=`wc -l raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt | cut -f 2 | sed 's/ .*//g'`; if (($l==0)); then echo "STARTING"; cat genolik.$chrom.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 63 $f.h 02.yellow.h $f | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; fi ; echo $chrom; done

## Add chrom and merge
for id in `seq 1 63`; do name=`head -$id ./refpanel.012.indv | tail -1`; for chrom in `seq 1 20`; do sed -i 's/^/\t/' ./raw_calls/full_refpanel.35kb.p2.$id.$chrom.txt; sed -i s/^/$chrom/ ./raw_calls/full_refpanel.35kb.p2.$id.$chrom.txt; done; cat ./raw_calls/full_refpanel.35kb.p2.$id.* > ./raw_calls/refpanel.full_refpanel.35kb.p2.$name.$cov.txt; done

mkdir tracts
module load R; for id in `seq 42 63`; do name=`head -$id ../refpanel.012.indv | tail -1`; sed -e s/INDIV/$name/g run.get_tracts.R > r.$id.R; sbatch r.$id.R; rm r.$id.R; done

