# Approximate the timing of admixture using DATES (from Priya Moorjani)
## Genotype and recombination files need to be updated to the most recent genome version

## We'll grab Amboseli genotype calls from: /data/tunglab/tpv/panubis1_genotypes/calls_merged/merged_shared.allchroms.vcf.gz
## Masked refpanel genotype calls from: /data/tunglab/tpv/panubis1_genotypes/masked_final/masked.filtered.no_chr.yes_intersect_50.vcf.gz
# High coverage yellow baboons list: 00_amboseli.list
# High coverage refpanel samples: 00_yel.list and 00_anu.list

cd /data/tunglab/tpv/local_ancestry/DATES
# Extract high coverage amboseli individuals
module load vcftools; module load plink
vcftools --gzvcf /data/tunglab/tpv/panubis1_genotypes/calls_merged/merged_shared.allchroms.vcf.gz --keep 00_amboseli.list --maf 0.01 --max-alleles 2 --recode --out hicov.amboseli
vcf-sort hicov.amboseli.recode.vcf > tmp.ambo.recode.vcf
## Manually remove underscores in sample names
# get rid of single nucleotide deletions, then convert to ancestrymap format
grep -v '^#' tmp.ambo.recode.vcf | cut -f 1-8 > tmp.ambo.snps; grep '\*' tmp.ambo.snps | cut -f 1-2 > tmp.ambo.to_exclude.list
vcftools --exclude-positions-overlap tmp.ambo.to_exclude.list --vcf tmp.ambo.recode.vcf --recode --out tmp.ambo.for_plink
plink --vcf tmp.ambo.for_plink.recode.vcf --snps-only --maf 0.05 --recode --out plink.ambo 
/data/tunglab/tpv/Programs/EIG-6.1.4/bin/convertf -p par.PED.ANCESTRYMAP2
## Manually edit the population column : ambo.v2.ind

# Extract high coverage, masked refpanel samples 
vcftools --gzvcf /data/tunglab/tpv/panubis1_genotypes/masked_final/masked.filtered.no_chr.yes_intersect_50.vcf.gz --keep 00_anu.list --keep 00_yel.list --recode --out hicov.ref
vcf-sort hicov.ref.recode.vcf > tmp.hicov.recode.vcf
## Manually remove underscores in sample names
# get rid of single nucleotide deletions, then convert to ancestrymap format
grep -v '^#' tmp.hicov.recode.vcf | cut -f 1-8 > tmp.hicov.snps; grep '\*' tmp.hicov.snps | cut -f 1-2 > tmp.ref.to_exclude.list
vcftools --exclude-positions-overlap tmp.ref.to_exclude.list --vcf tmp.hicov.recode.vcf --recode --out tmp.ref.for_plink
plink --vcf tmp.ref.for_plink.recode.vcf --snps-only --maf 0.05 --recode --out plink.n46 
/data/tunglab/tpv/Programs/EIG-6.1.4/bin/convertf -p par.PED.ANCESTRYMAP
## Manually edit the population column for yellow vs anubis : n46.v2.ind

## Reformat snp file, adding positions 
module load R; R
library(data.table); fread("ambo.snp") -> snp
#dim(snp); subset(snp, snp$V4 >= min(rcr$`#`) & snp$V4 <= max(rcr$left_snp)) -> snp ; dim(snp) #clear out SNPs beyond RCR rate
s2 <- as.data.frame(cbind(paste("snp", snp$V4,".", snp$V2, snp$V5, snp$V6,sep=""), snp$V2))
colnames(s2)[1:2] <- c("snpid","chr")
s2$GP <- NA
s2$PP <- snp$V4
load("/data/tunglab/tpv/Windows_25kb_recombination.RData")
subset(rcr, rcr$n24_anubis < 100*median(rcr$n24_anubis)) -> res2
s3 <- NULL 
for (i in unique(snp$V2)) {
	tmp=subset(s2, s2$chr == i) 
	r1=subset(res2, res2$chr == paste("chr",i,sep=""))
	r1$sum <- cumsum(r1$n24_anubis)
	## baboon chromsome lengths
	d <- c(172,115,100,164,113,111,116,93,88,88,91,126,77,76,69,89,77,84,63,63) [i]
	
	r1$cM <- r1$sum*d/max(r1$sum)
	lo <- loess(r1$cM ~ r1$end, span = 0.01)
	tmp$GP <- predict(lo, tmp$PP)
	tmp$GP[tmp$PP < 75000] <- predict(lo,75000)*tmp$PP[tmp$PP < 75000]/ 75000
	tmp <- tmp[!is.na(tmp$GP),]
	rbind(s3,tmp) -> s3; print(i)
}
#s2$GP <- s2$PP/1e6
write.table(s3, "ambo.v2.snp", row.names=F, col.names=F, sep="\t", quote=F)


library(data.table); fread("n46.snp") -> snp
s2 <- as.data.frame(cbind(paste("snp", snp$V4,".", snp$V2, snp$V5, snp$V6,sep=""), snp$V2))
colnames(s2)[1:2] <- c("snpid","chr")
s2$GP <- NA
s2$PP <- snp$V4
load("/data/tunglab/tpv/Windows_25kb_recombination.RData")
subset(rcr, rcr$n24_anubis < 100*median(rcr$n24_anubis)) -> res2
s3 <- NULL 
for (i in unique(snp$V2)) {
	tmp=subset(s2, s2$chr == i) 
	r1=subset(res2, res2$chr == paste("chr",i,sep=""))
	r1$sum <- cumsum(r1$n24_anubis)
	
	d <- c(172,115,100,164,113,111,116,93,88,88,91,126,77,76,69,89,77,84,63,63) [i]
	
	r1$cM <- r1$sum*d/max(r1$sum)
	lo <- loess(r1$cM ~ r1$end, span = 0.01)
	tmp$GP <- predict(lo, tmp$PP)
	tmp$GP[tmp$PP < 75000] <- predict(lo,75000)*tmp$PP[tmp$PP < 75000]/ 75000
	tmp <- tmp[!is.na(tmp$GP),]
	rbind(s3,tmp) -> s3
	print(i)
}
#s2$GP <- s2$PP/1e6
write.table(s3, "n46.v2.snp", row.names=F, col.names=F, sep="\t", quote=F)


## Merge ANCESTRYMAP files 
/data/tunglab/tpv/Programs/DATES-master/example/mergeit -p par.mergeit > mergeit.log 

## So there was an ordering issue with the imputation? 
library(data.table); fread("family_packed.snp") -> snp
fread("ambo.v2.snp") -> tst
tst$V3/100 -> tst$V3 # switch from cM to M 
subset(tst, tst$V1 %in% snp$V1) -> t2
t2$V3 -> snp$V3; rm(t2,tst)

snp$V1 <- paste("snp",row.names(snp),sep="")
s2 <- NULL 
for (i in unique(snp$V2)) {
	s=subset(snp, snp$V2 == i)
	s<-s[order(s$V4),]
	s$p <- c(-0.00025,s$V3[-nrow(s)])
	s$V3-s$p -> s$d 
	s$d[s$d < 1e-16] <- 1e-16
	s$V3 <- cumsum(s$d)
	#sum(order(s$V3)==order(s$V4))/nrow(s)
	rbind(s2,s) -> s2; print(i) 
}; write.table(s2[,1:6], "family_packed.snp", row.names=F, col.names=F, sep="\t", quote=F)

/data/tunglab/tpv/Programs/EIG-6.1.4/bin/convertf -p par.convertf.order


# Run DATES 
cd /data/tunglab/tpv/local_ancestry/DATES/
## Load relevant packages
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/tunglab/tpv/Programs/gsl-2.3/lib:/data/tunglab/tpv/Programs/fftw-3.3.3/lib
export LD_LIBRARY_PATH
export PATH=$PATH:$LD_LIBRARY_PATH
export PATH=$PATH:/data/tunglab/tpv/Programs/DATES-master/src/bin/
module load OpenBLAS/0.2.20-gcb01
module load glibc/2.14-gcb01
module load gnuplot 

/data/tunglab/tpv/Programs/DATES-master/src/bin/dates -p par.dates > log.dates 


