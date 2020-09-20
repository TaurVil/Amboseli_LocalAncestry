#### masked vs unmasked refpanel 
masking <- 'unmasked'

#########   Chunk to get data       ###########
### Recent vs historical
load(paste("~/Baboon/Resources/Recently_admixed_samples_",masking,".RData", sep=""))

a_data <- fread("~/Baboon/Paper1a_ancestrycalling/amboseli_tracts_fullrefpanel.txt")
colnames(a_data) <- c("chr", "start", "end", "state", "length", "u_prev", "u_after", "name")
a_data$length <- a_data$end - a_data$start

## Get "historical" tracts with between 10 and 500kb of homozygous anubis ancestry
t_hist_names <- ids$name[ids$historical == 'yes']
hist_tracts <- subset(a_data, a_data$name %in% t_hist_names & a_data$length > 10000 & a_data$length < 500000 & a_data$state == 2)
write.table(hist_tracts, "~/Baboon/Paper1a_ancestrycalling/historical_homozygous_ancestry.txt", row.names=F, col.names=T, sep="\t")

## Get "recent" tracts with >500kb of homozygous anubis ancestry
t_recent_names <- rownames(recent); t_recent_names <- ids$name[ids$sname %in% t_recent_names]
recent_tracts <- subset(a_data, a_data$name %in% t_recent_names & a_data$length > 500000 & a_data$state == 2)
write.table(recent_tracts, "~/Baboon/Paper1a_ancestrycalling/recent_homozygous_ancestry.txt", row.names=F, col.names=T, sep="\t")

## Each of these sum to ~19x individuals across the genome (assuming 3Gb which is a slight overestimate)
## Move to /data/tunglab/tpv/local_ancestry/
sed -i 's/\r//g' ../../local_ancestry/historical_hybrids.txt #Remove windows line breaks
sed -i 's/\r//g' ../../local_ancestry/recent_hybrids.txt

## Setup new working home
cd /data/tunglab/tpv/panubis1_genotypes/
mkdir amboseli_homozygous_segments; cd amboseli_homozygous_segments/

## Get bed file per individual of sites to maintain, and list of individuals for recent and historic
for f in `cat ../../local_ancestry/recent_hybrids.txt ../../local_ancestry/historical_hybrids.txt`; do grep $f ../../local_ancestry/amboseli_tracts_unmasked.txt > $f.homozygous.bed; done
## Remove half base pairs from the bed files
for f in `ls *bed`; do cat $f | cut -f 1-3 | sed 's/\.5//g' | sed 's/chr//g' > tmp; mv tmp $f; done
for f in `ls *bed`; do cat bed_head $f > tmp; mv tmp indiv_beds/$f; rm $f; done

mkdir indiv_vcfs
for f in `cat ../../local_ancestry/recent_hybrids.txt ../../local_ancestry/historical_hybrids.txt | sed 's/\r//g' `; do sed -e s/SAMPLE_NAME/$f/g run.01.get_indiv_vcf.sh > r.$f.sh; sbatch --mem=4G --nice r.$f.sh; echo $f; rm r.$f.sh; done


## Reform those individual files back into a full vcf, then filter for at least 10 individuals per call set
mkdir indiv_vcfs/recent; mkdir indiv_vcfs/historic
for f in `cat ../../local_ancestry/recent_hybrids.txt`; do mv indiv_vcfs/subvcf.$f*gz* indiv_vcfs/recent; done
for f in `cat ../../local_ancestry/historical_hybrids.txt`; do mv indiv_vcfs/subvcf.$f*gz* indiv_vcfs/historic; done

module load tabix; module load samtools; module load vcftools; module load bcftools
bcftools merge indiv_vcfs/historic/subvcf.*.recode.vcf.gz -O z -o ./historic_unfiltered.vcf.gz
bcftools merge indiv_vcfs/recent/subvcf.*.recode.vcf.gz -O z -o ./recent_unfiltered.vcf.gz

mkdir refpanel_vcfs; for vers in `cat 00_versions.list`; do sed -e s/VERSION_NAME/$vers/g run.02.merge_filter_vcfs.sh > r.$vers.sh; sbatch --mem=8G r.$vers.sh; rm r.$vers.sh; done 

