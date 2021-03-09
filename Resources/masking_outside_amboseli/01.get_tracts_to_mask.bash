## Goal: to mask variants which are likely heterospecific in baboons outside of Amboseli 

# Inputs: lclae_to_mask.txt and ibdmix_tracts.SPECIES.N.txt where species is yellow or anubis and N is 30, 50, or 70 and refers to the proportion of potential source individuals. These input files were saved locally as the consequence of XXX script. 

# make working directory on the cluster
cd /data/tunglab/tpv/local_ancestry/masking/ 
mkdir masking_tracts;

## Get IBDmix and LCLAE tracts 
## brought files over from local drive. Renamed lclae as lclae.txt from lclae_to_mask.txt 
## attach 'chr' start to the ibdmix tracts (which just have the chromosome number). LCLAE already has the `chr`
for f in `ls masking_tracts/ibdmix_tracts.*`; do sed -i s/^/chr/g $f; sed -i s/chrchr/chr/g $f; echo $f; done 

## Get subsets of lclae tracts for high coverage samples. This 00_anu.list and 00_yel.list refer to 28 and 18 individuals respectively, high coverage baboons from SNPRC founders (n=7 yellow; 24 anubis), the baboon genome project (anubis from Aberdare and SNPRC; 1 yellow from Mikumi), and Mikumi (n=10). 
## Get subsets for low coverage samples too. This is the individuals not in the previous set, and includes 4 Mikumi baboons, 7 Maasai Mara Baboons, and 6 WaNPRC baboons. 
for f in `cat /data/tunglab/tpv/panubis1_genotypes/high_coverage/00_anu.list`; do grep $f masking_tracts/lclae.txt >> masking_tracts/lclae_highcov.txt; echo $f; done 
for f in `cat /data/tunglab/tpv/panubis1_genotypes/high_coverage/00_yel.list`; do grep $f masking_tracts/lclae.txt >> masking_tracts/lclae_highcov.txt; echo $f; done 
cp masking_tracts/lclae.txt masking_tracts/lclae_lowcov.txt; for f in `cat /data/tunglab/tpv/panubis1_genotypes/high_coverage/00_anu.list`; do grep -v $f masking_tracts/lclae_lowcov.txt > tmp.txt; mv tmp.txt masking_tracts/lclae_lowcov.txt; echo $f; done; for f in `cat /data/tunglab/tpv/panubis1_genotypes/high_coverage/00_yel.list`; do grep -v $f masking_tracts/lclae_lowcov.txt > tmp.txt; mv tmp.txt masking_tracts/lclae_lowcov.txt; echo $f; done 

## use bedtools to get the filtering sets, which are combinations of LCLAE and IBDmix  
module load bedtools2; cut -f 8 masking_tracts/lclae.txt | sort | uniq | grep -v 'name' > 00_refpanel.names
### No low cov 
##### LCLAE 
bedtools sort -i masking_tracts/lclae_highcov.txt | sed -e 's/\.5//g'  > masking_tracts/no_lclae.bed 
##### IBDmix  
cat masking_tracts/ibdmix_tracts.*.30.txt | grep -v length > tmp.txt; bedtools sort -i tmp.txt > masking_tracts/no_ibd_30.bed ; rm tmp.txt 
cat masking_tracts/ibdmix_tracts.*.50.txt | grep -v length > tmp.txt; bedtools sort -i tmp.txt > masking_tracts/no_ibd_50.bed ; rm tmp.txt 
cat masking_tracts/ibdmix_tracts.*.70.txt | grep -v length > tmp.txt; bedtools sort -i tmp.txt > masking_tracts/no_ibd_70.bed ; rm tmp.txt 
##### union
for name in `cat 00_refpanel.names`; do cat masking_tracts/no_ibd_30.bed masking_tracts/no_lclae.bed | grep $name | sed -e 's/\.5//g' | cut -f 1-3 > tmp.txt; bedtools sort -i tmp.txt > tmp2.txt; bedtools merge -i tmp2.txt | sed "s/$/\t${name}/"  >> masking_tracts/none_union_30.bed; rm tmp.txt; rm tmp2.txt; echo $name; done 
for name in `cat 00_refpanel.names`; do cat masking_tracts/no_ibd_50.bed masking_tracts/no_lclae.bed | grep $name | sed -e 's/\.5//g' | cut -f 1-3 > tmp.txt; bedtools sort -i tmp.txt > tmp2.txt; bedtools merge -i tmp2.txt | sed "s/$/\t${name}/"  >> masking_tracts/none_union_50.bed; rm tmp.txt; rm tmp2.txt; echo $name; done 
for name in `cat 00_refpanel.names`; do cat masking_tracts/no_ibd_70.bed masking_tracts/no_lclae.bed | grep $name | sed -e 's/\.5//g' | cut -f 1-3 > tmp.txt; bedtools sort -i tmp.txt > tmp2.txt; bedtools merge -i tmp2.txt | sed "s/$/\t${name}/"  >> masking_tracts/none_union_70.bed; rm tmp.txt; rm tmp2.txt; echo $name; done 
##### intersect
for name in `cat 00_refpanel.names`; do cat masking_tracts/ibdmix_tracts.*.30.txt | grep $name > tmp.ibd.bed; grep $name masking_tracts/lclae_highcov.txt | sed -e 's/\.5//g' > tmp.lclae.bed; bedtools intersect -a tmp.ibd.bed -b tmp.lclae.bed | sed "s/$/\t${name}/"  >> masking_tracts/none_intersect_30.bed ; rm tmp.ibd.bed; rm tmp.lclae.bed; echo $name; done 
for name in `cat 00_refpanel.names`; do cat masking_tracts/ibdmix_tracts.*.50.txt | grep $name > tmp.ibd.bed; grep $name masking_tracts/lclae_highcov.txt | sed -e 's/\.5//g' > tmp.lclae.bed; bedtools intersect -a tmp.ibd.bed -b tmp.lclae.bed | sed "s/$/\t${name}/"  >> masking_tracts/none_intersect_50.bed ; rm tmp.ibd.bed; rm tmp.lclae.bed; echo $name; done 
for name in `cat 00_refpanel.names`; do cat masking_tracts/ibdmix_tracts.*.70.txt | grep $name > tmp.ibd.bed; grep $name masking_tracts/lclae_highcov.txt | sed -e 's/\.5//g' > tmp.lclae.bed; bedtools intersect -a tmp.ibd.bed -b tmp.lclae.bed | sed "s/$/\t${name}/"  >> masking_tracts/none_intersect_70.bed ; rm tmp.ibd.bed; rm tmp.lclae.bed; echo $name; done 

### with low cov 
##### LCLAE 
cat masking_tracts/lclae.txt | sed -e 's/\.5//g' > masking_tracts/yes_lclae.bed 
##### IBDmix  
cat masking_tracts/lclae_lowcov.txt masking_tracts/no_ibd_30.bed | sed -e 's/\.5//g' > masking_tracts/yes_ibd_30.bed 
cat masking_tracts/lclae_lowcov.txt masking_tracts/no_ibd_50.bed | sed -e 's/\.5//g' > masking_tracts/yes_ibd_50.bed 
cat masking_tracts/lclae_lowcov.txt masking_tracts/no_ibd_70.bed | sed -e 's/\.5//g' > masking_tracts/yes_ibd_70.bed 
##### union
cat masking_tracts/none_union_30.bed masking_tracts/lclae_lowcov.txt | sed -e 's/\.5//g' > masking_tracts/yes_union_30.bed 
cat masking_tracts/none_union_50.bed masking_tracts/lclae_lowcov.txt | sed -e 's/\.5//g' > masking_tracts/yes_union_50.bed 
cat masking_tracts/none_union_70.bed masking_tracts/lclae_lowcov.txt | sed -e 's/\.5//g' > masking_tracts/yes_union_70.bed 
##### intersect
cat masking_tracts/none_intersect_30.bed masking_tracts/lclae_lowcov.txt | sed -e 's/\.5//g' > masking_tracts/yes_intersect_30.bed 
cat masking_tracts/none_intersect_50.bed masking_tracts/lclae_lowcov.txt | sed -e 's/\.5//g' > masking_tracts/yes_intersect_50.bed 
cat masking_tracts/none_intersect_70.bed masking_tracts/lclae_lowcov.txt | sed -e 's/\.5//g' > masking_tracts/yes_intersect_70.bed 

## get list of version to run 
ls masking_tracts/ | grep 'bed' | sed 's/.bed//g' > 00_versions.list


## At this point, let's clean up intermediary files. 
### we'll remove lclae*.txt and ibdmix_tracts.*txt from masking_tracts/, so we're left with just the bed files 
