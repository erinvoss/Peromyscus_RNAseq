#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --partition=savio2
#SBATCH --account=fc_nachman
#SBATCH --array=1-60
#SBATCH --output=slurmout/fastqc_test_job_%A_task_%a.out
#SBATCH --error=slurmout/fastqc_test_job_%A_task_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erinvoss@berkeley.edu

cd /global/scratch/users/erinvoss/Peromyscus-RNA

module load fastqc

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
inpath="01-FastP_Preproc"
outpath="analysis"

echo "SAMPLE: ${sample}"

fastqc ${inpath}/${sample}/${sample}_R1_cleaned.fastq.gz ${inpath}/${sample}/${sample}_R2_cleaned.fastq.gz \
-o ${outpath}/fastqc_reports_2022_10_17 -t 10

end=`date +%s`
runtime=$((end-start))
echo $runtime