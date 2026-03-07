#!/bin/bash

#SBATCH --job-name=sqlite-download
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=savio
#SBATCH --account=fc_nachman
#SBATCH --output=slurmout/trinotate_download_job_%A.out
#SBATCH --error=slurmout/trinotate_download_job_%A.err

cd /global/scratch/users/erinvoss/Peromyscus-RNA/03-Trinotate-Annotation

module load python/3.9.12
source activate trinity_env

TRINOTATE_HOME="/global/scratch/users/erinvoss/.conda/envs/trinity_env_2/bin"
transcriptome_inpath="/global/scratch/users/erinvoss/Peromyscus-RNA/02-Trinity-Assembly/PECA-Trinity-Output"
database_inpath="/global/scratch/users/erinvoss/Peromyscus-RNA/03-Trinotate-Annotation/Trinotate_Databases"
resultsDir="/global/scratch/users/erinvoss/Peromyscus-RNA/03-Trinotate-Annotation/PECA-Annotation"
shortName="Trinity-cd-hit-95" 

start=`date +%s`
echo $HOSTNAME

# Initiate database
echo "Initiating $shortName database"
Trinotate --db PECA_95_Trinotate.sqlite --init \
           --gene_trans_map ${transcriptome_inpath}/$shortName".gene_trans_map" \
           --transcript_fasta ${transcriptome_inpath}/$shortName".fa" \
           --transdecoder_pep ${database_inpath}/$shortName".transdecoder.pep"

# Load sequence analysis results into SQLite database
echo "Loading Sequence Analyses into Trinotate Database"
Trinotate PECA_95_Trinotate.sqlite LOAD_swissprot_blastx ${database_inpath}/$shortName".blastx.outfmt6" # Load transcript hits 
Trinotate PECA_95_Trinotate.sqlite LOAD_swissprot_blastp ${database_inpath}/$shortName".outfmt6" # Load protein hits
Trinotate PECA_95_Trinotate.sqlite LOAD_custom_blast --outfmt6 ${database_inpath}/$shortName"_custom_pema_ensembl_blastx.outfmt6" \
    --prog blastx --dbtype custom_PEMA_genome_blastx # Load custom P. maniculatus ensembl blast hits
Trinotate PECA_95_Trinotate.sqlite LOAD_custom_blast --outfmt6  ${database_inpath}/$shortName"_custom_peca_ncbi_blastx.outfmt6" --prog blastx \
	--dbtype custom_PECA_genome # Load custom P. californicus ncbi blast hits
Trinotate PECA_95_Trinotate.sqlite LOAD_pfam ${database_inpath}/$shortName"peca_95_pfam.out" # Load load Pfam domain entries
Trinotate PECA_95_Trinotate.sqlite LOAD_signalp ${database_inpath}/$shortName"peca_95_signalp.out" # Load signal peptide predictions
Trinotate PECA_95_Trinotate.sqlite LOAD_tmhmm ${database_inpath}/$shortName"tmhmm.out" # Load transmembrane domains

# Run Trinotate main step: generate annotation report 
# Optional: Use low E-value to filter stringently for high-quality hits
echo "Generating Trinotate Annotation Report"
Trinotate PECA_95_Trinotate.sqlite report -E 10e-10 --incl_trans > ${resultsDir}/$shortName"_annotation_report.xls"

# Generate annotation statistics report
echo "Generating Summary Statistics from Trinotate Annotation Report"
${TRINOTATE_HOME}/count_table_fields.pl ${resultsDir}/$shortName"_annotation_report.xls" > ${resultsDir}/$shortName"_annotation_report_stats.xls"

# Make Gene-Transcript Map 
$TRINITY_HOME/support_scripts/get_Trinity_gene_to_trans_map.pl \
    ${data_inpath}/$shortName".fasta.transdecoder.cds.fasta" > $shortName".fasta.gene_trans_map"        

# Make Annotation Feature Map using custom P. maniculatus ENSEMBL hits
# ERV Note: I worked through all the sequence analysis tools recommended in the Trinotate GitHub / Wiki, but after a detailed comparison of outputs, 
# I decided to identify transcripts according to the P. maniculatus ENSEMBL CDS database as the best balance of completeness and reproducibility. 
echo "Making annotation feature map using custom P. maniculatus ENSEMBL hits"
${TRINOTATE_HOME}/Trinotate_get_feature_name_encoding_attributes_custom_PEMA.pl \ 
                    # Note: this script was modified to map annotations to P. maniculatus ensembl identifiers
                 ${resultsDir}/$shortName"_annotation_report.xls"> ${resultsDir}/$shortName"_annot_feature_map.txt"             


end=`date +%s`
runtime=$((end-start))
echo $runtime