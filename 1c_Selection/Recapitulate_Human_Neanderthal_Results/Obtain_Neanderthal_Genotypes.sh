## Download Neanderthal genotypes from Svante Paabo's group

for f in `seq 1 22`; do wget http://ftp.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.${f}.mod.vcf.gz; echo $f; done
for f in `seq 1 22`; do wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr${f}.noRB.vcf.gz; mv chr${f}.noRB.vcf.gz Chagyrskaya.${f}.vcf.gz; echo $f; done
for f in `seq 1 22`; do wget http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr${f}_mq25_mapab100.vcf.gz; mv chr${f}_mq25_mapab100.vcf.gz Vindija.${f}.vcf.gz

module load vcftools; module load tabix
for f in `seq 1 22`; do gunzip Chagyrskaya.${f}.vcf.gz ; bgzip Chagyrskaya.${f}.vcf; tabix Chagyrskaya.${f}.vcf.gz; bcftools annotate -x ^FORMAT/GT Chagyrskaya.${f}.vcf.gz -I -O z -o Chagyrskaya.${f}.small.vcf.gz; done
for f in `seq 1 22`; do tabix Vindija.${f}.vcf.gz; bcftools annotate -x ^FORMAT/GT Vindija.${f}.vcf.gz -I -O z -o Vindija.${f}.small.vcf.gz; done

zcat AltaiNea.hg19_1000g.${f}.mod.vcf.gz | grep -v 'LowQual' | grep -v '\./\.' | bgzip > Altai.${f}.no_lowqual.vcf.gz; tabix Altai.${f}.no_lowqual.vcf.gz; bcftools annotate -x INFO,^FORMAT/GT Altai.${f}.no_lowqual.vcf.gz -I -O z –o Altai.${f}.small.vcf.gz

zcat AltaiNea.hg19_1000g.${f}.mod.vcf.gz | grep -v 'LowQual' | bgzip > Altai.${f}.no_lowqual.vcf.gz
tabix Altai.${f}.no_lowqual.vcf.gz
bcftools annotate -x INFO,^FORMAT/GT Altai.${f}.no_lowqual.vcf.gz -I -O z –o

tabix Chagyrskaya.$f.vcf.gz
bcftools annotate -x INFO,^FORMAT/GT Chagyrskaya.$f.vcf.gz -I -O v -o Chagyrskaya.small.$f.vcf
bgzip Chagyrskaya.small.$f.vcf; tabix Chagyrskaya.small.$f.vcf.gz 

f=22; bcftools merge Altai.$f.short.vcf.gz Chagyrskaya.small.$f.vcf.gz Vindija.$f.small.vcf.gz -O v -o Neanderthal.$f.vcf; vcftools --vcf Neanderthal.$f.vcf --freq --max-missing=0.5 --out N.$f; sed -i 's/:/\t/g' N.$f.frq
