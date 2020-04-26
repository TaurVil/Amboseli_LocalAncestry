# LocalAncestry
Pipelines for calling local ancestry from simulated or resequencing data. 

### Part 1: Map reads using bowtie2 to the new Wall assembly. 

##Get SRA files: 
for f in `cat 00_Wall2016_SRA`; do sed -e s/SRRNUMBER/$f/g /data/tunglab/tpv/scripts/download_SRR_SRA_fastq.sh > g.sh; sbatch g.sh; done

sed -e s/SRRNUMBER/SRR2565914/g ../scripts/download_SRR_SRA_fastq.sh > g.sh; sbatch g.sh
sed -e s/SRRNUMBER/SRR2650000/g ../scripts/download_SRR_SRA_fastq.sh > g.sh; sbatch g.sh
sed -e s/SRRNUMBER/SRR2565912/g ../scripts/download_SRR_SRA_fastq.sh > g.sh; sbatch g.sh
sed -e s/SRRNUMBER/SRR2650075/g ../scripts/download_SRR_SRA_fastq.sh > g.sh; sbatch g.sh


#Map the 290 low coverage sequencing files Arielle produced
cd /data/tunglab/tpv/LocalAncestry/ ; for f in `cat Amboseli_290/00samples2`; do sed -e s/FILE/$f/g ./run.01.mapBowtie2_wall.FILE.sh > g.sh; sbatch --mem=26000 --nice g.sh; done; rm g.sh

# bowtie2 -p 16 -t -x ./bowtie2_Wall15944/bowtie2_Wall15944 -1 /data/tunglab/asf40/wgs_data/MedGenome_ftp/FASTQs/FILE_R1_*.fastq.gz -2 /data/tunglab/asf40/wgs_data/MedGenome_ftp/FASTQs/FILE_R2_*.fastq.gz -S mapped.FILE.wall.sam --no-unal
# samtools view -bS mapped.FILE.wall.sam > Amboseli_290/mapped.FILE.wall.bam
# rm mapped.FILE.wall.sam



#Saved output files to /data/tunglab/tpv/LocalAncestry/Amboseli_290/output_mapping/
#Saves mapped bam files to /data/tunglab/tpv/LocalAncestry/Amboseli_290/, under the name of mapped.NAME.wall.bam

### Part 2: sort and add RG
cd /data/tunglab/tpv/LocalAncestry/Amboseli_290/ ; for f in `cat 00samples2`; do sed -e s/NUMBER/$f/g ../run.02.sort.FILE.sh > g.sh; sbatch --mem=16000 g.sh; done; rm g.sh

for f in `cat 00samples2`; do sed -e s/SAMPLE/$f/g ../run.03.addRG.FILE.sh > g.sh; sbatch --mem=16000 --nice g.sh; done; rm g.sh

#call genotypes. 
sbatch --array=1-274%100 --mem=26000 run.01.GATK_wall.sh
#for f in `cat 00scaffolds`; do for h in `cat 00covs`; do cat run.small_GATK.sh | sed -e s/CHROM/$f/g | sed -e s/COVS/$h/g > s.$f.$h.sh; sbatch --mem=34000 s.$f.$h.sh; done ; done 
#for f in `cat 00scaffolds`; do for h in `cat 00covs`; do cat run_merge_vcfs.sh | sed -e s/CHROM/$f/g | sed -e s/COVS/$h/g > s.$f.$h.sh; sbatch --mem=34000 s.$f.$h.sh; done ; done 



for f in `cat 00_wall_BigScaffolds`; do mv Amboseli290.$f.filt3.$f.vcf.gz ready.$f.vcf.gz ; done 
ls ready.Super-Scaffold_*.vcf.gz > 00ready
sed -i s/ready.//g 00ready
sed -i s/.vcf.gz//g 00ready


##Get genolik matrix
for f in `cat 00ready`; do cat get_genolik.sh | sed -e s/SCAFFOLD/$f/g  > s.$f.sh; sbatch --mem=8000 s.$f.sh; chmod 777 s.$f.sh; done 


##Get calls for each individual in a window (g) 
for f in `cat 00ready`; do for g in `seq 1 5`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 6 15`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 16 35`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 36 55`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 56 95`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 96 135`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 136 155`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 156 175`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 176 205`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 206 235`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 236 255`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 256 285`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 286 305`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done
for f in `cat 00ready`; do for g in `seq 306 344`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done

##Fix those which may not have finished. Iterate until complete (slurm output lengths are 0)
###Change precurser number for values 1 to 9
##started at 4. Done with 4, 5, 6, 7 , 8, 9, 3, 
wc -l 9*.*35kb*.txt > n_calls
sed -i 's/^[ \t]*//' n_calls
grep '^0' n_calls > no_calls
sed -i 's/\./\t/g' no_calls
wc -l no_calls 

num=`wc -l no_calls | awk '{print $1}'`
for i in `seq 1 $num`; do sed "${i}q;d" no_calls > tmp; f=`awk '{print $2}' tmp`; g=`awk '{print $5}' tmp`;  sed -e s/SCAFFOLD/$g/g get_callsh.sh | sed -e s/NUMBER/$f/g > g.$f.$g.sh; sbatch --out=out.$f.$g.out g.$f.$g.sh; rm g.$f.$g.sh; done

##Let's find a way to do this based on failed ones from the output file... 
wc -l out.*Super*.out > n_out
sed -i 's/^[ \t]*//' n_out 
grep '^3' n_out > failed 
sed -i 's/\./\t/g' failed
wc -l failed
rm out.*.out
n=`wc -l failed | awk '{print $1}'`
for i in `seq 1 $n`;  do sed "${i}q;d" failed > tmp2; f=`awk '{print $3}' tmp2`; g=`awk '{print $4}' tmp2`;  sed -e s/SCAFFOLD/$g/g get_callsh.sh | sed -e s/NUMBER/$f/g > g.$f.$g.sh; sbatch --out=out.$f.$g.out --mem=5500 g.$f.$g.sh; rm g.$f.$g.sh; done


#Get LCLAE calls & information 
for f in `cat /data/tunglab/tpv/LocalAncestry/ALLBAMS/00_wall_BigScaffolds`; do cat get_genolik.sh | sed -e s/SCAFFOLD/$f/g  > s.sh; sbatch --mem=26000 s.sh; done ; rm s.sh 
	for f in `cat 00ready`; do cat get_genolik.sh | sed -e s/SCAFFOLD/$f/g  > s.$f.sh; sbatch --mem=8000 s.$f.sh; chmod 777 s.$f.sh; done ; 


for f in `cat 00ready`; do for g in `seq 101 125`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done

for f in `cat 00ready`; do for g in `seq 293 344`; do cat get_callsh.sh | sed -e s/SCAFFOLD/$f/g | sed -e s/NUMBER/$g/g > g.sh; sbatch --out out.$f.$g.out g.sh; rm g.sh; done ; done


wc -l *35kb*.txt > n_calls
sed -i 's/^[ \t]*//' n_calls
grep '^0' n_calls > no_calls
sed -i 's/\./\t/g' no_calls
wc -l no_calls 

for i in `seq 1 1500`; do sed "${i}q;d" no_calls > tmp; f=`awk '{print $2}' tmp`; g=`awk '{print $5}' tmp`;  sed -e s/SCAFFOLD/$g/g get_callsh.sh | sed -e s/NUMBER/$f/g > g.$f.$g.sh; sbatch g.$f.$g.sh; rm g.$f.$g.sh; done


#add chromosome name to each file
 for h in `seq 1 344`; do for g in `cat 00ready`; do sed -i 's/^/\t/' $h.35kb.d2.$g.txt ; sed -i s/^/$g/g $h.35kb.d2.$g.txt; done; echo $h; done
 
 for h in `seq 1 5`; do cat $h.35kb.d2.Super* > $h.35kb.d2.txt; done


for f in `seq 1 344`; do cat run_mode_min50.R | sed -e s/INDIV/$f/g > g.$f.sh; sbatch --mem=45000 --nice g.$f.sh; rm g.$f.sh; done; 
for f in `seq 293 344`; do cat run_mode_min50.R | sed -e s/INDIV/$f/g > g.$f.sh; sbatch --mem=45000 --nice g.$f.sh; rm g.$f.sh; done; 


















