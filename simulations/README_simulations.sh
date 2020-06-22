#README_simulations.txt
# SELAM instruction manual: https://github.com/russcd/SELAM/blob/master/SELAM_manual.pdf
cd /data/tunglab/tpv/local_ancestry/simulated_data/

###Run SELAM
#Load in modules and pre-requisites
module load gcc; module load samtools; module load python/2.7.6-fasrc01; module load virtualenv; module load gsl; 
# First time we will need to create the virtual environment: virtualenv venv/
source /data/tunglab/tpv/local_ancestry/simulated_data/venv/bin/activate; pip install --upgrade pip==6.0.7
# and install SELAM: wget https://github.com/russcd/SELAM/archive/master.zip; unzip master.zip; cd /data/tunglab/tpv/LocalAncestry/SELAM_simulated_data/SELAM-master/src; make

#So let's output 4 chromosomes for each individual. 
/data/tunglab/tpv/Programs/SELAM/src/SELAM -d try1_demography.txt -o try1_output.txt --seed 112 -c 5 77 84 63 63 1
#-c says to call 4 chromosomes, with the lengths given in morgans, which correspond to the smallest chromosomes in the baboon genome (17-20).
# The last chromosome is inherited from only the maternal line, hence why it has a weird output. 
deactivate

ls i*output.txt | sed 's/.output.txt//g' > 00names 

## Let's clean up the tracks a little bit, removing the comment lines and the last chromosome (#4), which is from the maternal lineage only
for f in `ls i*output.txt`; do grep -v '^#' $f | grep -v -P "\t4\t" > tmp; mv tmp $f; done

##We must convert ancestry output into tracts for the homozygous individuals
mkdir tracts; module load R; for f in `cat 00names `; do sed -e s/NAME/$f/g run.get_tracts.R > s.$f.sh; sbatch s.$f.sh; rm s.$f.sh; done
## There are checks in place to make sure each call is biallelic and length > 0

##Generate a vcf file for each sample
# creat example header to use
zcat /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed/02.yel.20.recode.vcf.gz | grep '^#' | tail -1 > vcf_example_header
# manually edit to create 1 name only ("SAMPLE")

# get yellow and anubis allele frequencies
module load vcftools; cd /data/tunglab/tpv/panubis1_genotypes/calls_unadmixed; for f in `seq 17 20`; do vcftools --gzvcf 02.yel.$f.recode.vcf.gz --freq --out chr$f.yellow; sed -i '1d' chr$f.yellow.frq; sed -i 's/:/\t/g' chr$f.yellow.frq; vcftools --gzvcf 02.anu.$f.recode.vcf.gz --freq --out chr$f.anubis; sed -i '1d' chr$f.anubis.frq; sed -i 's/:/\t/g' chr$f.anubis.frq; done
mv chr*frq /data/tunglab/tpv/local_ancestry/simulated_data/; cd /data/tunglab/tpv/local_ancestry/simulated_data/; 

# get a simplified vcf file for each sample
mkdir simulated_vcfs; module load R; for f in `cat 00names`; do for g in `cat 00chroms`; do cat run.get_sample_vcf.R | sed -e s/NAME/$f/g | sed -e s/SCAF/$g/g > g.$f.$g.R; sbatch g.$f.$g.R; rm g.$f.$g.R; done; done 

# We need to add the full vcf header to each sample
zcat /data/tunglab/tpv/panubis1_genotypes/calls_merged/04.merged_shared.1.vcf.gz | grep '^#' | sed '$ d' > full_header.vcf
for f in `cat 00names`; do cat full_header.vcf > simulated_vcfs/$f.vcf; head -1 simulated_vcfs/$f.chr20.vcf >> simulated_vcfs/$f.vcf; for g in `ls simulated_vcfs/$f.chr*vcf`; do sed '1d' $g >> simulated_vcfs/$f.vcf; done; done

# Get small genome to generate reads for
module load samtools; for f in `cat 00chroms`; do samtools faidx /data/tunglab/shared/genomes/panubis1/Panubis1.0.fa $f > $f.fa; echo $f; done; cat chr*.fa > Reduced_Genome.fa; rm chr*.fa

##Simulate reads using NEAT-genreads, a tool developed by Zachary Stephens 
#github: https://github.com/zstephens/neat-genreads
#publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5125660/
# Generate 10x coverage in SE, 100bp reads
mkdir sim_reads
for f in `cat 00names`; do sed s/INDIV/$f/g run.get_fasta.sh > g.$f.sh; sbatch g.$f.sh; done





##We need to modify these a bit to get only SNPs of interest and include a proper header  
for f in `ls i*Super-Scaff*.vcf`; do sed '1d' $f > 2.$f; done 
for f in `ls 2.i*`; do sed -e s/FILE/$f/g get_simple_vcf.R > g.$f.sh; sbatch g.$f.sh; done

rm 2.i*.vcf ; mv i*Super*.vcf True_VCFs/
for f in `ls 3.2.*`; do cat header.vcf $f > temp; mv temp $f; echo $f; done 


