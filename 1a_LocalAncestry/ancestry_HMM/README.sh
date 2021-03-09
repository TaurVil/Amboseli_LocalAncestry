## Run Ancestry HMM for all Amboseli individuals

cd /data/tunglab/tpv/local_ancestry/HMM

### Prep data from vcf file 
module load bcftools; bcftools query -f '%CHROM \t %POS \t %REF \t %ALT [\t %AD]\n' /data/tunglab/tpv/panubis1_genotypes/masked_final/final.amboseli_with_masked_refpanel.vcf.gz > temp
sed -i 's/\./0,0/g' temp
sed -i 's/,/\t/g' temp

module load bcftools ; bcftools query -l /data/tunglab/tpv/panubis1_genotypes/masked_final/final.amboseli_with_masked_refpanel.vcf.gz > samples.list

cut -f 1-2 temp > temp2.sites 
## There are 445 Amboseli individuals, which are at the front of the file. They will therefore be columsn 5:894 in the data 
cut -f 5-894 temp > temp2.ambodata
cut -f 895-1020 temp > temp2.refdata; grep -v 'AMB' samples.list > refsamples.list 
## There are 63 reference panel samples 
## yellow samples are 2, 6, 7, 22, 26, 27, 28, 39-52, 57: so that's 3:4, 11:14, 43:44, 51:56, 77:104, 113:114
## We only want the SW yellow which are in the first 31, so that's 3:4, 11:14, 43:44, 51:56
## Rest are anubis, so we want those in the first 62 columns

## so in R, let's get the sum of anubis ref, anubis alt, yellow ref, and yellow alt allele counts. 
module load R; R
library(data.table); 
fread("temp2.sites", header=F) -> sites
fread("temp2.refdata", header=F) -> data
colnames(sites) <- c('chrom', 'snp')

yel <- data[,c(3:4, 11:14, 43:44, 51:56, 77:104, 113:114)]
anu <- data[,!c(3:4, 11:14, 43:44, 51:56, 77:104, 113:114)]

keep <- seq(1,48,2); sites$ref_anu <- rowSums(anu[,..keep], na.rm=T)
keep <- seq(2,48,2); sites$alt_anu <- rowSums(anu[,..keep], na.rm=T)
keep <- seq(1,14,2); sites$ref_yel <- rowSums(yel[,..keep], na.rm=T)
keep <- seq(2,14,2); sites$alt_yel <- rowSums(yel[,..keep], na.rm=T)
## calculate distance between SNPs assuming a uniform recombination rate of 1 cM/Mb (or 1 M/1e8 bp). Express in M/bp. 
sites$dist <- sites$snp - c(1, sites$snp[-nrow(sites)]) 
sites$dist <- sites$dist * 1e-8

## export intermediary without amboseli data. This was done because we were exceeding memory trying to do it all at once. 
write.table(, "data_no_amboseli.txt", row.names=F, col.names=T, sep="\t", quote=F)

## add in ambodata
fread("temp2.ambodata") -> ambo; 
cbind(sites, ambo) -> sites
write.table(sites, "data.panel", row.names=F, col.names=F, sep="\t", quote=F)

## Prepare to run Ancestry HMM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/tunglab/tpv/Programs/conda/lib/
export PATH=$PATH:/data/tunglab/tpv/Programs/conda/lib/
module load Anaconda; conda --version
module load libtool; libtool --version

# Install the AncestryHMM software
# wget https://github.com/russcd/Ancestry_HMM/archive/master.zip; gunzip master.zip
# cd Ancestry_HMM/src/
# make



## Run Ancestry HMM for each chromosome
for chrom in `cat 00chroms`; do mkdir $chrom; cp test.samples $chrom/ ; grep '^${chrom}' data.panel > $chrom/data.panel ; done
## for chrom 1 and 2, manually subset the rows of data.panel to avoid having other chromosomes which start with the same number 
## remove sites that are not differentiated between populations
module load R; R
library(data.table); fread("data.panel") -> d
yel <- d$V4/rowSums(d[,3:4])
anu <- d$V6/rowSums(d[,5:6])
d2 <- subset(d, abs(anu-yel) > 0.2)
write.table(d2, "filtered.panel", row.names=F, col.names=F, sep="\t", quote=F)

/data/tunglab/tpv/Programs/Ancestry_HMM/Ancestry_HMM-master/src/ancestry_hmm -i filtered.panel -s test.samples -a 2 0.25 0.75 -p 1 100000 0.75 -p 0 -100 0.25 --fixed



# Create ancestry tracts from posterior calls
# array_for_get_ancestry.sh with INDIV and CHROM, which calls get_ancestryHMM_tracts.R 
for chrom in `cat 00chroms`; do sed -e s/CHROM/$chrom/g get_anceestryHMM_tracts.R > get.$chrom.R; sed -e s/CHROM/$chrom/g array_for_get_ancestry.sh > array.$chrom.sh; sbatch --array=1-445 array.$chrom.sh; rm array.$chrom.sh; done
# outputs a file for CHROM_INDIV in tracts

# Convert tracts to HMM tracts
for i in `cat 00_ambo_names.list`; do cat tracts/*${i}* > HMMtracts/${i}.txt; done 

## We'll pull these to the local directory to analyze
