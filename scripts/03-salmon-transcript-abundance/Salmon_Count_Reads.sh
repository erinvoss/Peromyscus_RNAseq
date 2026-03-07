#!/bin/bash

#SBATCH --job-name=salmon-peca
#SBATCH --nodes=1
#SBATCH --time=36:00:00
#SBATCH --partition=savio2_htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=fc_nachman
#SBATCH --output=slurmout/salmon-abundance-peca_job_%A.out
#SBATCH --error=slurmout/salmon-abundance-peca_job_%A.err

start=`date +%s`
echo $HOSTNAME

cd /global/scratch/users/erinvoss/Peromyscus-RNA/04B-Salmon-Abundance/Salmon-Spring-2023/salmon-enspem-transcriptome-annot
fastqDir='/global/scratch/users/erinvoss/Peromyscus-RNA/01-FastP_Preproc'
transcriptomeDir='/global/scratch/users/erinvoss/Peromyscus-RNA/03-Trinotate-Annotation/Annot_Transcriptomes_02_21_2023'
salmonDir='/global/scratch/users/erinvoss/Peromyscus-RNA/04B-Salmon-Abundance/Salmon-Spring-2023/salmon-enspem-transcriptome-annot'
TRINITY_HOME='/global/home/users/erinvoss/Modules/sources/trinityrnaseq-v2.15.1'

module load python/3.9.12
source activate salmon_env

# Main Salmon Transcript Abundance step
echo "Estimating Transcript Abundance with Salmon"
while read file in files; do
	echo "Processing $file%{.txt}"
    ${TRINITY_HOME}/util/align_and_estimate_abundance.pl \ # Trinity built in script to call salmon to quantify transcript abundance 
                   --transcripts ${transcriptomeDir}/PECA_transcriptome_ENSPEM_annot_final.fa \ # Transcriptome with transcripts renamed with Ensembl identifiers
                   --prep_reference \
                   --seqType fq \ # raw RNAseq data is in fastq format
                   --samples_file ${salmonDir}/$file \ # Specifies location of fastq files for each sample to count reads 
                   --est_method salmon \ # Use salmon to quantify
                   --gene_trans_map ${transcriptomeDir}/PEBO_transcriptome_ENSPEM_annot_final.gene_trans_map \ # Use this file to pair transcripts with genes
                   --output_dir ${salmonDir}/PEBO --thread_count 4;
done < PECA_sample_lists.txt 
	# Samples grouped according to species and tissue; all replicates of one tissue processed together in Salmon
	# Repeat this step with each of three species for which ERV has RNAseq data and a transcriptome 
	# After Salmon, need to rename output files (quant.sf) according to species, tissue, and replicate to prep for import into R and DeSeq2. 
     

runtime=$((end-start))
echo $runtime


