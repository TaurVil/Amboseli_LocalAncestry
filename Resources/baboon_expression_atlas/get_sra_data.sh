#!/bin/bash
#SBATCH --get-user-env

### sbatch --array=1-767 get_sra_data.sh


module load samtools
module load sratoolkit
module load tabix 

index=${SLURM_ARRAY_TASK_ID}
srr=`head -$index 00_SRRsamples.list | tail -1`



prefetch $srr
vdb-validate $srr
fasterq-dump $srr

bgzip $srr.fastq







