# Make merged refpanel of yellow and anubis samples, for all chromosomes

module load bcftools

for f in `seq 1 20`; do bcftools merge /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed/02.anu.$f.recode.vcf.gz /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed/02.yel.$f.recode.vcf.gz -O z -o /data/tunglab/tpv/local_ancestry/$f.vcf.gz; echo $f; done

bcftools concat /data/tunglab/tpv/local_ancestry/*.vcf.gz -O z -o /data/tunglab/tpv/local_ancestry/refpanel.vcf.gz

vcftools --gzvcf ../refpanel.vcf.gz --mac 2 --max-missing 0.5 --kept-sites --out ../refpanel

## Add "chr" to kept-sites file & vcf
sed -e 's/^\([0-9XY]\)/chr\1/' ../refpanel.kept.sites > ../refpanel.kept_chroms.sites; mv ../refpanel.kept_chroms.sites ../refpanel.kept.sites

zcat ../refpanel.vcf.gz | sed -e 's/##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' -e 's/^\([0-9XY]\)/chr\1/' > ../tmp.vcf
mv ../tmp.vcf ../refpanel.vcf; bgzip ../refpanel.vcf; tabix ../refpanel.vcf.gz
