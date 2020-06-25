## Can we detect introgression in our reference panel samples? 

# Get genotype likelihoods for unadmixed individuals
## The chromosome level vcfs are housed in /data/tunglab/tpv/local_ancestry/
cd /data/tunglab/tpv/local_ancestry/; module load vcftools; vcftools --gzvcf 20.vcf.gz --012 --out refpanel

mkdir unadmixed_individuals; cd unadmixed_individuals
sbatch --array=1-20 --mem=5G run.get_genolik.sh

## 41 anubis 1-41
### SW: 1-24, Mara 25-31, Tulane 32-37, BGP 38-41
## yellow 42-63
### SW: 42-48, Mikumi 49-63

mkdir raw_calls
echo 41 `seq 1 41` > 02.anubis.h; echo 22 `seq 42 63` > 02.yellow.h
echo 21 `seq 42 63` > 02.yellow_minus1.h

sbatch --array=44-62 run.get_raw_calls.sh

## 63 is special and needs to be done manually

for f in `echo 42`; do sed -e "s/ $f / /g" 02.yellow_minus1.h > $f.h; name=`head -$f ../refpanel.012.indv | tail -1`; for chrom in `seq 1 20`; do touch ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; l=`wc -l raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt | cut -f 2 | sed 's/ .*//g'`; if (($l==0)); then echo "STARTING"; cat genolik.$chrom.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 64 02.anubis.h $f.h $f | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; fi ; echo $chrom; done; echo $f; rm $f.h; done


for f in `seq 42 63`; do sed -e "s/ $f / /g" 02.yellow_minus1.h > $f.h; name=`head -$f ../refpanel.012.indv | tail -1`; for chrom in `seq 1 20`; do cat genolik.$chrom.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 64 02.anubis.h $f.h $f | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; echo $chrom; done; echo $f; rm $f.h; done

for f in `echo 43`; do sed -e "s/ $f / /g" 02.yellow_minus1.h > $f.h; name=`head -$f ../refpanel.012.indv | tail -1`; for chrom in `seq 1 20`; do cat genolik.$chrom.genolik | /data/tunglab/tpv/Programs/LCLAE/filtbaboon2c 64 02.anubis.h $f.h $f | /data/tunglab/tpv/Programs/LCLAE/geno_lik2 .2 35000 > ./raw_calls/full_refpanel.35kb.p2.$f.$chrom.txt; echo $chrom; done; echo $f; rm $f.h; done

