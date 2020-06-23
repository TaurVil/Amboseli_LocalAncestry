# Make merged refpanel of yellow and anubis samples, for all chromosomes

module load bcftools

for f in `seq 1 20`; do bcftools merge /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed/02.anu.$f.recode.vcf.gz /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed/02.yel.$f.recode.vcf.gz -O z -o /data/tunglab/tpv/local_ancestry/$f.vcf.gz; echo $f; done

bcftools concat /data/tunglab/tpv/local_ancestry/*.vcf.gz -O z -o /data/tunglab/tpv/local_ancestry/refpanel.vcf.gz

