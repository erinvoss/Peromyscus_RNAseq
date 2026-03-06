#!/bin/bash

#SBATCH --job-name=fastp
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --partition=savio2
#SBATCH --account=fc_nachman
#SBATCH --array=1-60
#SBATCH --output=slurmout/fastp_job_%A_task_%a.out 
#SBATCH --error=slurmout/fastp_job_%A_task_%a.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erinvoss@berkeley.edu

cd /global/scratch/users/erinvoss/Peromyscus-RNA

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

inpath="00-RawData"
outpath="01-FastP_Preproc"
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load fastp

call="fastp \
  -i ${inpath}/${sample}/${sample}_R1.fastq.gz \
  -I ${inpath}/${sample}/${sample}_R2.fastq.gz \
  -o ${outpath}/${sample}/${sample}_R1_cleaned.fastq.gz \
  -O ${outpath}/${sample}/${sample}_R2_cleaned.fastq.gz \
  -n 5 \
  -q 20 \
  -u 30 \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_window_size=4 \
  --cut_mean_quality=20 \
  --length_required=25 \
  -j ${outpath}/${sample}/${sample}_fastp_report.json \
  -h ${outpath}/${sample}/${sample}_fastp_report.html \
  -w 16"


echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
