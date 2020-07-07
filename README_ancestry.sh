# OUTDATED BITS
Pipelines for calling local ancestry from simulated or resequencing data. 

#call genotypes. 
sbatch --array=1-274%100 --mem=26000 run.01.GATK_wall.sh
#for f in `cat 00scaffolds`; do for h in `cat 00covs`; do cat run.small_GATK.sh | sed -e s/CHROM/$f/g | sed -e s/COVS/$h/g > s.$f.$h.sh; sbatch --mem=34000 s.$f.$h.sh; done ; done 
#for f in `cat 00scaffolds`; do for h in `cat 00covs`; do cat run_merge_vcfs.sh | sed -e s/CHROM/$f/g | sed -e s/COVS/$h/g > s.$f.$h.sh; sbatch --mem=34000 s.$f.$h.sh; done ; done 

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


