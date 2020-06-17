#!/usr/bin/env Rscript
#SBATCH --get-user-env

library(data.table); library(plyr)
fread("./tracts/NAME.tracts.txt") -> sample 
read.delim("vcf_example_header") -> VCF

subset(sample, sample$chrom == "SCAF") -> sample
	
nom=paste("SCAF",".anubis.frq",sep="")
anu <- fread(nom); rm(nom); anu$V1 <- paste("chr",anu$V1,sep="")
colnames(anu) <- c('chr', 'pos','n_alleles','n_chr','ref','ref_freq','alt','alt_freq')
nom=paste("SCAF",".yellow.frq",sep="")
yel <- fread(nom); rm(nom); yel$V1 <- paste("chr",yel$V1,sep="")
colnames(yel) <- c('chr', 'pos','n_alleles','n_chr','ref','ref_freq','alt','alt_freq')

vcf <- as.data.frame(matrix(ncol=10, nrow=nrow(anu))); colnames(vcf) <- colnames(VCF)
k=1; unique(c(anu$pos, yel$pos)) -> sites

set.seed(100)

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

				if (temp$ancestry[1] == 0) {
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
				if (temp$ancestry[1] == 2) {
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
				if (temp$ancestry[1] == 1) {
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

write.table(vcf, "simulated_vcfs/NAME.SCAF.vcf", row.names=F, col.names=T, quote=F, sep="\t")


