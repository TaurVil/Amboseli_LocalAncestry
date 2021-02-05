# Approximate the timing of admixture using DATES (from Priya Moorjani)
## Genotype and recombination files need to be updated to the most recent genome version

# Run DATES 
## Load relevant packages
LD_LIBRARY_PATH=/data/tunglab/tpv/Programs/gsl-2.3/lib:/data/tunglab/tpv/Programs/fftw-3.3.3/lib
export LD_LIBRARY_PATH
export PATH=$PATH:$LD_LIBRARY_PATH
export PATH=$PATH:/data/tunglab/tpv/dating_admixture/DATES-master/src/bin
module load OpenBLAS/0.2.20-gcb01
module load glibc/2.14-gcb01
module load gnuplot 


## We'll grab Amboseli genotype calls from: /data/tunglab/tpv/panubis1_genotypes/calls_merged/merged_shared.allchroms.vcf.gz
## Masked refpanel genotype calls from: /data/tunglab/tpv/panubis1_genotypes/masked_final/masked.filtered.no_chr.yes_intersect_50.vcf.gz
# High coverage yellow baboons list: 00_amboseli.list
# High coverage refpanel samples: 00_yel.list and 00_anu.list

cd /data/tunglab/tpv/local_ancestry/DATES
# Extract high coverage amboseli individuals
module load vcftools; module load plink
vcftools --gzvcf /data/tunglab/tpv/panubis1_genotypes/calls_merged/merged_shared.allchroms.vcf.gz --keep 00_amboseli.list --maf 0.01 --max-alleles 2 --recode --out hicov.amboseli
plink --vcf hicov.amboseli.recode.vcf --maf 0.05 --recode --out plink.ambo 
/data/tunglab/tpv/Programs/EIG-6.1.4/bin/convertf -p par.PED.ANCESTRYMAP2
## Manually edit the population column 

# Extract high coverage, masked refpanel samples 
vcftools --gzvcf /data/tunglab/tpv/panubis1_genotypes/masked_final/masked.filtered.no_chr.yes_intersect_50.vcf.gz --keep 00_anu.list --keep 00_yel.list --recode --out hicov.ref
plink --vcf hicov.ref.recode.vcf --maf 0.05 --recode --out plink.n33
/data/tunglab/tpv/Programs/EIG-6.1.4/bin/convertf -p par.PED.ANCESTRYMAP
## Manually edit the population column for yellow vs anubis 

## Reformat snp file, adding positions 
/data/tunglab/tpv/SW_anubis_founders/n24_anubis_output


library(data.table); fread("ambo.snp") -> snp
#dim(snp); subset(snp, snp$V4 >= min(rcr$`#`) & snp$V4 <= max(rcr$left_snp)) -> snp ; dim(snp) #clear out SNPs beyond RCR rate
s2 <- as.data.frame(cbind(paste("snp", snp$V4,sep=""), snp$V2))
colnames(s2)[1:2] <- c("snpid","chr")
s2$GP <- NA
s2$PP <- snp$V4
load("~/n10_recombination_rates.RData")
subset(res, res$n24_anubis < 100*median(res$n24_anubis)) -> res2
s3 <- NULL 
for (i in unique(snp$V2)) {
	tmp=subset(s2, s2$chr == i) 
	r1=subset(res2, res2$chr == paste("chr",i,sep=""))
	r1$sum <- cumsum(r1$n24_anubis)
	
	if (i == 1) {d <- 280}; if (i == 2) {d <- 264}; if (i == 3) {d <- 220}; if (i == 4) {d <- 210}; if (i ==5) {d <- 200}
	if (i == 6) {d <- 190}; if (i == 7) {d <- 185}; if (i == 8) {d <- 170}; if (i == 9) {d <- 165}; if (i ==10) {d <- 175}
	if (i == 11) {d <- 160}; if (i == 12) {d <- 175}; if (i == 13) {d <- 130}; if (i == 14) {d <- 125}; if (i ==15) {d <- 130}
	if (i == 16) {d <- 130}; if (i == 17) {d <- 130}; if (i == 18) {d <- 120}; if (i == 19) {d <- 108}; if (i ==20) {d <- 108}
	
	r1$cM <- r1$sum*d/max(r1$sum)
	lo <- loess(r1$cM ~ r1$end, span = 0.01)
	tmp$GP <- predict(lo, tmp$PP)
	tmp$GP[tmp$PP < 250000] <- predict(lo,250000)*tmp$PP[tmp$PP < 250000]/ 250000
	rbind(s3,tmp) -> s3
}
#s2$GP <- s2$PP/1e6
write.table(s3, "ambo.v2.snp", row.names=F, col.names=F, sep="\t", quote=F)


library(data.table); fread("n33.snp") -> snp
#dim(snp); subset(snp, snp$V4 >= min(rcr$`#`) & snp$V4 <= max(rcr$left_snp)) -> snp ; dim(snp) #clear out SNPs beyond RCR rate
s2 <- as.data.frame(cbind(paste("snp", snp$V4,sep=""), snp$V2))
colnames(s2)[1:2] <- c("snpid","chr")
s2$GP <- NA
s2$PP <- snp$V4
load("~/n10_recombination_rates.RData")
subset(res, res$n24_anubis < 100*median(res$n24_anubis)) -> res2
s3 <- NULL 
for (i in unique(snp$V2)) {
	tmp=subset(s2, s2$chr == i) 
	r1=subset(res2, res2$chr == paste("chr",i,sep=""))
	r1$sum <- cumsum(r1$n24_anubis)
	
	if (i == 1) {d <- 280}; if (i == 2) {d <- 264}; if (i == 3) {d <- 220}; if (i == 4) {d <- 210}; if (i ==5) {d <- 200}
	if (i == 6) {d <- 190}; if (i == 7) {d <- 185}; if (i == 8) {d <- 170}; if (i == 9) {d <- 165}; if (i ==10) {d <- 175}
	if (i == 11) {d <- 160}; if (i == 12) {d <- 175}; if (i == 13) {d <- 130}; if (i == 14) {d <- 125}; if (i ==15) {d <- 130}
	if (i == 16) {d <- 130}; if (i == 17) {d <- 130}; if (i == 18) {d <- 120}; if (i == 19) {d <- 108}; if (i ==20) {d <- 108}
	
	r1$cM <- r1$sum*d/max(r1$sum)
	lo <- loess(r1$cM ~ r1$end, span = 0.01)
	tmp$GP <- predict(lo, tmp$PP)
	tmp$GP[tmp$PP < 250000] <- predict(lo,250000)*tmp$PP[tmp$PP < 250000]/ 250000
	#tmp$GP[tmp$PP > max(r1$end)] <- predict(lo,250000)*tmp$PP[tmp$PP < 250000] d/250000
	## Just did a little repair instead for chrom 8 where this failed: subset out that chromosome, took the mean rate across sites (pbp), assigned for  the sites with NAs (s$V3[is.na(s$V3)] <- pbp*s$V4[is.na(s$V3)]), then insert back into the full matrix (snp$V3[snp$V2 == 8] <- s$V3). save output again 
	print("binding")
	rbind(s3,tmp) -> s3
	print(i)
}
#s2$GP <- s2$PP/1e6
write.table(s3, "n33.v2.snp", row.names=F, col.names=F, sep="\t", quote=F)


## Merge ANCESTRYMAP files 
/data/tunglab/tpv/dating_admixture/DATES-master/example/mergeit -p par.mergeit > mergeit.log 

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


## Run Dates 
/data/tunglab/tpv/dating_admixture/DATES-master/src/bin/dates -p par.dates > log.dates 


