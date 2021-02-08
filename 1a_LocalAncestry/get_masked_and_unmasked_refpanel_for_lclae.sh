### Merge masked and unmasked genotypes for the reference panel samples. 
## masked versions have the "masked_" prefix added to each name 
## final file is masked_and_unmasked_refpanel_samples/masked_and_unmasked.vcf.gz 

cd /data/tunglab/tpv/panubis1_genotypes
mkdir masked_and_unmasked_refpanel_samples

cp masked_final/masked.filtered.yes_intersect_50.vcf.gz masked_and_unmasked_refpanel_samples/masked.vcf.gz 
gunzip masked_and_unmasked_refpanel_samples/masked.vcf.gz 

cd masked_and_unmasked_refpanel_samples/
grep '^#' masked.vcf > head.vcf
grep -v '^#' masked.vcf > rest.vcf; sed -i 's/^chr//' rest.vcf 

vi masked_and_unmasked_refpanel_samples/masked.vcf ## add "masked" to the name of each individual

cat head.vcf  rest.vcf > masked_named.vcf
bgzip masked_named.vcf 

rm masked.vcf; rm head.vcf; rm rest.vcf 

cd ..
tabix masked_and_unmasked_refpanel_samples/masked_named.vcf.gz
module load bcftools; bcftools merge  masked_and_unmasked_refpanel_samples/masked_named.vcf.gz anubis.vcf.gz yellow.vcf.gz -O z -o masked_and_unmasked_refpanel_samples/masked_and_unmasked.vcf.gz 


cat calls_merged/03.shared.*sites > masked_and_unmasked_refpanel_samples/sites_to_keep

vcftools --gzvcf masked_and_unmasked.vcf.gz  --positions sites_to_keep --recode --out masked_and_unmasked_variable

