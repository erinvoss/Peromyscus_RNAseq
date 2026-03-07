#!/bin/bash

#SBATCH --job-name=cd-hit-est
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --partition=savio
#SBATCH --account=fc_nachman
#SBATCH --output=slurmout/peca_transcriptome_post_job_%A.out
#SBATCH --error=slurmout/peca_transcriptome_post_job_%A.err

cd /global/scratch/users/erinvoss/Peromyscus-RNA/02-Trinity-Assembly/PECA-Trinity-Output

module load python/3.9.12
source activate trinity_env

# Step 1: collapse transcripts with greater than 95% identity together 

echo 'Removing redundant transcripts with cd-hit-est'

cd-hit-est -i Trinity.fasta -o Trinity-cd-hit-95 -c 0.95 -p 1 -d 0 -b 3 -M 64000 -T 0

# Step 2: filter out any chimeras, mis-assemblies etc with RNAquast, compared with PECA reference genome

cd /global/scratch/users/erinvoss/Peromyscus-RNA/02-Trinity-Assembly/PECA-Trinity-Output/PECA-Trinity-95

RNAQUAST_PATH=/global/home/users/erinvoss/Modules/sources/rnaQUAST-2.2.2

referenceDir='/global/scratch/users/erinvoss/Peromyscus-RNA/References/Genomes/P_californicus'
transcriptomeDir='/global/scratch/users/erinvoss/Peromyscus-RNA/02-Trinity-Assembly/PECA-Trinity-Output/PECA-Trinity-95'

echo 'Assessing PECA Transcriptome Quality with rnaQuast'

python ${RNAQUAST_PATH}/rnaQUAST.py \
 	--transcripts ${transcriptomeDir}/Trinity-cd-hit-95.fa  \
    --reference ${referenceDir}/GCF_007827085.1_ASM782708v3_genomic.fna \
    --gtf ${referenceDir}/GCF_007827085.1_ASM782708v3_genomic.gtf \
    --output_dir ${outputDir}/peca_ncbi_rnaQUAST_results

