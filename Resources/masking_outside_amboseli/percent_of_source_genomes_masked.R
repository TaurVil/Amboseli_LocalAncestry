## Getting the proportion of each reference panel genome that was masked (SW founder individuals)

## get just the hi coverage samples for removed tracts 
## head -209643 masking_tracts/yes_intersect_50.bed > hi_cov.tracts.removed.bed

library(data.table); fread("hi_cov.tracts.removed.bed") -> data
anu <- fread("/data/tunglab/tpv/unadmixed_individuals/00_anu.list", header=F)[1:24,]
yel <- fread("/data/tunglab/tpv/unadmixed_individuals/00_yel.list", header=F)[11:17,]
chroms <- fread("/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa.fai")

total_length <- sum(chroms$V2[1:20])
anu$masked <- 0; for (i in 1:nrow(anu)) {
  tmp <- subset(data, data$V6 == anu$V1[i])
  tmp$V5 <- tmp$V3 - tmp$V2
  anu$masked[i] <- sum(tmp$V5)
}
yel$masked <- 0; for (i in 1:nrow(yel)) {
  tmp <- subset(data, data$V6 == yel$V1[i])
  tmp$V5 <- tmp$V3 - tmp$V2
  yel$masked[i] <- sum(tmp$V5)
}

mean(anu$masked)/total_length; sd(anu$masked)/total_length
min(anu$masked)/total_length; max(anu$masked)/total_length

mean(yel$masked)/total_length; sd(yel$masked)/total_length
min(yel$masked)/total_length; max(yel$masked)/total_length
