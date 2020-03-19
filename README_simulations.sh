#README_simulations.txt

###Run SELAM
#Load in modules and pre-requisites
module load gcc; module load samtools; module load python/2.7.6-fasrc01; module load virtualenv; module load gsl; 
source /data/tunglab/tpv/LocalAncestry/SELAM_simulated_data/venv/bin/activate; pip install --upgrade pip==6.0.7

#So let's output 9 chromosomes for each individual. The chromosome lengths correspond to scaffolds in the Wall assembly, such that the lengths of ancestry tracts in the simulated data resemble the expected distribution in real data. 
/data/tunglab/tpv/LocalAncestry/SELAM_simulated_data/SELAM-master/src/SELAM -d try1_demography.txt -o try1_output.txt --seed 112 -c 9 3 7 7 11 4 10 8 15 12 
#-c says to call 9 chromosomes, with the lengths given in morgans. so we'll hopefully be getting a lot of recombination events.

##We must convert  to the ancestry tract lengths we get 
module load R 
for f in `cat 00names `; do sed -e s/NAME/$f/g get_tracts.R > s.$f.sh; sbatch s.$f.sh; done
##This script currently codes each 


##Generate a vcf file for each sample  
#Read in tracts
#Read in anubis/yellow frequencies 
####for f in `cat 00names`; do cat get_vcfs.R | sed -e s/NAME/$f/g  > g.$f.R; sbatch g.$f.R; done 
for f in `cat 00names`; do for g in `cat 00_selam_chroms`; do cat get_single_vcf.R | sed -e s/NAME/$f/g | sed -e s/SCAF/$g/g > g.$f.$g.R; sbatch g.$f.$g.R; done; done 

#Outputs vcf file for each individual 
##We need to modify these a bit to get only SNPs of interest and include a proper header  
for f in `ls i*Super-Scaff*.vcf`; do sed '1d' $f > 2.$f; done 
for f in `ls 2.i*`; do sed -e s/FILE/$f/g get_simple_vcf.R > g.$f.sh; sbatch g.$f.sh; done

rm 2.i*.vcf ; mv i*Super*.vcf True_VCFs/
for f in `ls 3.2.*`; do cat header.vcf $f > temp; mv temp $f; echo $f; done 


