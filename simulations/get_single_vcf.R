#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table); library(plyr)
fread("NAME_tracts.txt") -> sample 
read.delim("vcf_example_header") -> VCF

fread("freq_SCAF_anu.frq") -> data
colnames(data) <- c('chr', 'pos','n_alleles','n_chr','ref','alt')

vcf <- as.data.frame(matrix(ncol=10, nrow=nrow(data)))
colnames(vcf) <- colnames(VCF)
k=1

strsplit(data$ref,':') -> r1
strsplit(data$alt,':') -> a1
ldply(r1) -> df1
ldply(a1) -> df2
cbind(data, df1, df2) -> data 
colnames(data) <- c('chr', 'pos', 'n_alleles', 'n_chr', 'r', 'a', 'ref', 'ref_freq', 'alt', 'alt_freq')
write.table(data, "t1.NAME.SCAF", row.names=F, col.names=T, sep="\t", quote=F)
fread("t1.NAME.SCAF") -> anu
	
fread("freq_SCAF_yel.frq") -> data
colnames(data) <- c('chr', 'pos','n_alleles','n_chr','ref','alt')
strsplit(data$ref,':') -> r1
strsplit(data$alt,':') -> a1
ldply(r1) -> df1
ldply(a1) -> df2
cbind(data, df1, df2) -> data 
colnames(data) <- c('chr', 'pos', 'n_alleles', 'n_chr', 'r', 'a', 'ref', 'ref_freq', 'alt', 'alt_freq')
write.table(data, "t1.NAME.SCAF", row.names=F, col.names=T, sep="\t", quote=F)
fread("t1.NAME.SCAF") -> yel
	
rm(data) 
	
	unique(c(anu$pos, yel$pos)) -> sites 


subset(sample, sample$chrom == "SCAF") -> sample 
	
for (j in 1:length(sites)) {
	#For each site, first get the ancestry status
	subset(sample, sample$start <= sites[j] & sample$end > sites[j]) -> temp 
	subset(anu, anu$pos == sites[j] ) -> tanu
	subset(yel, yel$pos == sites[j]) -> tyel
	if (nrow(tanu) == 1 & nrow(tyel) == 1) {
		if (tanu$ref == tyel$ref & tanu$alt == tyel$alt & !(is.na(tyel$ref_freq[1])) & !(is.na(tanu$ref_freq[1])) ) {
			#Populate the SNP characteristics
				vcf$CHROM[k] <- tanu$chr[1]
				vcf$POS[k] <- tanu$pos[1]
				vcf$REF[k] <- tanu$ref [1]
				vcf$ALT[k] <- tanu$alt [1]
				
				if (temp$status[1] == 0) {
					#make two draws from the yellow baboons
					if (tyel$ref_freq == 1) { vcf$sampleINDIV[k] <- "0/0"} 
					if (tyel$alt_freq == 1) {vcf$sampleINDIV[k] <- "1/1"}
					if (tyel$ref_freq < 1 & tyel$alt_freq < 1) {
						runif(2,0,1) -> d; sum(d > tyel$ref_freq) -> jk
						if (jk == 0) {vcf$sampleINDIV[k] <- "0/0"}
						if (jk == 1) {vcf$sampleINDIV[k] <- "0/1"}
						if (jk == 2) {vcf$sampleINDIV[k] <- "1/1"}
					}
				}
				if (temp$status[1] == 2) {
					#make two draws from the anubis baboons
					if (tanu$ref_freq == 1) { vcf$sampleINDIV[k] <- "0/0"} 
					if (tanu$alt_freq == 1) {vcf$sampleINDIV[k] <- "1/1"}
					if (tanu$ref_freq < 1 & tanu$alt_freq < 1) {
						runif(2,0,1) -> d; sum(d > tanu$ref_freq[1]) -> jk
						if (jk == 0) {vcf$sampleINDIV[k] <- "0/0"}
						if (jk == 1) {vcf$sampleINDIV[k] <- "0/1"}
						if (jk == 2) {vcf$sampleINDIV[k] <- "1/1"}
					}
				}
				if (temp$status[1] == 1) {
					#make one draw from each species 
					runif(1,0,1) -> d; sum(d > tanu$ref_freq[1]) -> a
					runif(1,0,1) -> d; sum(d > tyel$ref_freq[1]) -> y
					sum(a,y) -> jk 
						if (jk == 0) {vcf$sampleINDIV[k] <- "0/0"}
						if (jk == 1) {vcf$sampleINDIV[k] <- "0/1"}
						if (jk == 2) {vcf$sampleINDIV[k] <- "1/1"}
					}
				k <- k+1; print(k)
		}
	}
}
				

vcf$INFO <- vcf$ID <- "."
vcf$QUAL <- 100
vcf$FILTER <- "QD"
vcf$FORMAT <- "GT"

write.table(vcf, "NAME.SCAF.vcf", row.names=F, col.names=T, quote=F, sep="\t")


