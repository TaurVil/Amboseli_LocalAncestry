library(data.table); library(plyr)

## Identify sites that are not variable within modern African populations
african <- fread("1k_only_afr.txt")
a <- subset(african, !african$V9 %in% c("0", "0,0", "0,0,0") & african$V1 %in% 1:22)

a_quin <- a[grepl(a$V9,pattern=",.*,.*,"),]
a <- a[!grepl(a$V9,pattern=",.*,.*,"),]
a_quad <- a[grepl(a$V9,pattern=",.*,"),]
a <- a[!grepl(a$V9,pattern=",.*,"),]
a_tri <- a[grepl(a$V9,pattern=","),]
a <- a[!grepl(a$V9,pattern=","),]
strsplit(a_quin$V9, ",") -> V1; ldply(V1) -> V1; colnames(V1) <- c("alt1", "alt2", "alt3", "alt4"); cbind(a_quin, V1) -> a_quin; a_quin$nvar <- 4
strsplit(a_quad$V9, ",") -> V1; ldply(V1) -> V1; colnames(V1) <- c("alt1", "alt2", "alt3"); cbind(a_quad, V1) -> a_quad; a_quad$nvar <- 3; a_quad$alt4 <- NA
strsplit(a_tri$V9, ",") -> V1; ldply(V1) -> V1; colnames(V1) <- c("alt1", "alt2"); cbind(a_tri, V1) -> a_tri; a_tri$nvar <- 2; a_tri$alt4 <- a_tri$alt3 <- NA
a$V9 -> a$alt1; a$nvar <- 1; a$alt4 <- a$alt3 <- a$alt2 <- NA

rbind(a, a_tri, a_quad, a_quin) -> afr 
afr$min <- apply(afr[,10:14],1,min, na.rm=T)
afr$min <- as.numeric(as.character(afr$min))
dim(afr); afr2 <- subset(afr, afr$min > 0); dim(afr2)
write.table(afr2, "1k_africa_variable.txt", row.names=F, col.names=F, sep="\t", quote=F)

## bash: get head of neanderthal frequency data
## for f in `seq 1 22`; do cut -f 1-6 N.$f.frq > tmp.$f.frq; echo $f; done 

## switch back to R, for Neanderthal data
## cd /data/tunglab/tpv/my_genomes/humans/gr37/Neanderthal_vcfs/
## module load R; R
library(data.table); library(plyr)
fread("1k_africa_variable.txt") -> african

all_fixed <- NULL; for (k in 1:22) {
d <- fread(paste('tmp.',k,'.frq', sep=""))
colnames(d)[6] <- "ref_freq"
subset(d, d$ref_freq == 0) -> fixed 
a <- subset(african, african$V1 == k)

tmp <- subset(fixed, ! fixed$POS %in% a$V4)
rbind(all_fixed, tmp) -> all_fixed; rm(tmp)
}

for f in `seq 1 22`; do sed '1d' N.$f.frq | cut -f 2 > tmp.$f.pos; awk 'BEGIN{
     for(r=getline; r>0;){
         for(s=e=$1; (r=getline)>0 && $1<=e+1; e=$1);
         print s==e ? s : s","e
     }
     exit -r
 }' tmp.$f.pos > tmp.$f.intervals
 grep -v ',' tmp.$f.intervals | awk '{print $1, $1}' | sed 's/ /\t/g' > tmp2; grep ',' tmp.$f.intervals | sed -e 's/,/\t/g' > tmp1; cat tmp1 tmp2 | sort > tmp.$f.intervals; done
 
## in R we need to add one to the second column to get a length from it 


### R fragment. Not sure what the entire input is?

df4 <- d %>%
mutate( parts = strsplit(tmp,",")) %>%
group_by( V9 ) %>%
   do(  data.frame(
     {
       idx <- 1:length(.$parts[[1]])
       lst <- lapply(idx,
                     function(x) .$parts[[1]][x])
       names(lst) <- lapply(idx,
                            function(x) paste("allele",x,sep="") )

       (lst)
     } , stringsAsFactors=FALSE)
   ) %>%
 inner_join(d,by="V9")

df4$A1 <- as.numeric(as.character(df4$allele1))
df4$A2 <- as.numeric(as.character(df4$allele2))

