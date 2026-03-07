#!/bin/bash

#SBATCH --job-name=peca-busco
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=savio2
#SBATCH --account=fc_nachman
#SBATCH --output=slurmout/busco_peca_job_%A.out
#SBATCH --error=slurmout/busco_peca_job_%A.err

cd $PEROMYSCUS_HOME/02-Trinity-Assembly/PECA-Trinity-Output/PECA-Trinity-95

source activate busco_env

start=`date +%s`
echo $HOSTNAME

echo "Assessing P. californicus transcriptome completeness with BUSCHO Euarchontoglires"
# Repeat with additional transcriptomes as needed

busco -i Trinity-cd-hit-95.fa \
      -l /global/scratch/users/erinvoss/Peromyscus-RNA/References/Busco-Resources/euarchontoglires_odb10 \
      -o PECA_Busco_Output -m transcriptome --cpu 24 --offline

echo "Generating summary statistics for Trinity output transcriptome"
python /global/scratch/users/erinvoss/Peromyscus-RNA/scripts/trinity_scripts/trinity_contig_stats.pl Trinity-cd-hit-95.fa

end=`date +%s`
runtime=$((end-start))
echo $runtime
