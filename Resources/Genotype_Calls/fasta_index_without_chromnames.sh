## index no chromosome name file

cd /data/tunglab/tpv/panubis1_genotypes/

sed 's/chr//g' /data/tunglab/shared/genomes/panubis1/Panubis1.0.fa > ./Panubis1.nochromname.fa
java -jar ~/picard.jar CreateSequenceDictionary R=./Panubis1.nochromname.fa O=Panubis1.nochromname.dict

module load samtools; samtools faidx ./Panubis1.nochromname.fa

