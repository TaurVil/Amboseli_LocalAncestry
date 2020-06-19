## Working with 7 high coverage amboseli individuals, the bams for which are in /data/tunglab/tpv/mapped_bams/amboseli_hicov/

cd /data/tunglab/tpv/local_ancestry/downsampling


## Get coverage for each sample
ls /data/tunglab/tpv/mapped_bams/amboseli_hicov/ | grep -v 'bai' | sed s/_MarkDuplicates.bam//g > 00_names.txt 
for f in `cat 00_names.txt`; do sed -e s/FILENAME/$f/g run.calculate_coverage.sh > tmp.sh; sbatch --mem=8G --out cov.$f.out tmp.sh; rm tmp.sh; done

## Downsample to get 10x data, as well as 5x, 2x, 1x, and p5x. 
cat *out | sed 's/Average =  //g' > 00_covs.txt
sbatch --array=1-7 --mem=8G run.downsample.sh
# check the output for coverage information
for f in `cat 00_bams.txt`; do sed -e s/FILENAME/$f/g run.calculate_coverage.sh > tmp.sh; sbatch --mem=8G --out cov.$f.out tmp.sh; rm tmp.sh; done

## Convert bams to chrXX rather than just XX
for file in AMB*.bam; do filename=`echo $file | sed s/.bam//g`; samtools view -H $file | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > ${filename}_chr.bam; echo $file; echo ${filename}_chr.bam; done

## Call genotypes for each set of bams
mkdir gVCF; sbatch --array=1-35 --mem=16G  run.gvcf.sh
sbatch --array=1-7 --mem=16G  run.gvcf_original.sh
