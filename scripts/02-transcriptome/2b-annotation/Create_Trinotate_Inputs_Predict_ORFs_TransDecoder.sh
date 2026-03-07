#!/bin/bash

#SBATCH --job-name=transdecoder-peca
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --partition=savio
#SBATCH --account=fc_nachman
#SBATCH --output=slurmout/transdecoder_predict_job_%A.out
#SBATCH --error=slurmout/transdecoder_predict_job_%A.err


cd /global/scratch/users/erinvoss/Peromyscus-RNA/03-Trinotate-Annotation/PECA-Annotation/Trinotate-cd-hit-95

data_inpath="/global/scratch/users/erinvoss/Peromyscus-RNA/02-Trinity-Assembly/PECA-Trinity-Output"
transdecoder_dir="/global/scratch/users/erinvoss/Peromyscus-RNA/03-Trinotate-Annotation/PECA-Annotation/Trinotate-cd-hit-95/Trinity-cd-hit-95.fa.transdecoder_dir"
db_inpath="/global/scratch/users/erinvoss/Peromyscus-RNA/03-Trinotate-Annotation/Trinotate-Databases"
shortName="Trinity-cd-hit-95" 


start=`date +%s`

module load python/3.9.12
source activate trinity_env_2

# Comment or uncomment each sequence analysis tool as needed below to build Trinotate input data. 

# Build Trinotate SQLite Database

#TRINOTATE_HOME="/global/scratch/users/erinvoss/.conda/envs/trinity_env_2/bin"
#${TRINOTATE_HOME}/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate

#gunzip Pfam-A.hmm.gz
#hmmpress Pfam-A.hmm

# Extract ORFs > 100 amino acids long from transcriptome
# At this point, option to filter for "best" or longest ORF
TransDecoder.LongOrfs -t ${data_inpath}/$shortName".fa"

# Use BLASTP homology search to extract additional proteins with Uniprot database
blastp -query ${transdecoder_dir}/$shortName".transdecoder.pep" \
       -db ${db_inpath}/uniprot_sprot.pep  -max_target_seqs 1 \
       -outfmt 6 -evalue 1e-5 -num_threads 20 > $shortName".blastp.outfmt6"

# Use BLASTX homology search to extract additional transcripts with Uniprot database
blastx -query ${transdecoder_dir}/$shortName".fa" \
       	-db ${db_inpath}/uniprot_sprot.pep  -max_target_seqs 1 \
       	-outfmt 6 -evalue 1e-5 -num_threads 24 > $shortName".blastx.outfmt6"


# Use BLASTX homology search to extract additional transcripts with custom Peromyscus maniculatus peptide database
blastx -query ${data_inpath}/$shortName".fa" \
		-db ${db_inpath}/Genomes/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.pep.all.fa \
      	-num_threads 24 -max_target_seqs 1 -evalue 1e-5 \
      	-outfmt 6 > $shortName"_custom_pema_ensembl_blastx.outfmt6"

# Use BLASTX homology search to extract additional transcripts with custom Peromyscus californicus peptide database
blastx -query ${data_inpath}/$shortName".fa"  \
       -db ${db_inpath}/Genomes/GCF_007827085.1_ASM782708v3_protein.faa \
       -num_threads 24 -max_target_seqs 1 -evalue 1e-5 \
       -outfmt 6 > $shortName"_custom_peca_ncbi_blastx.outfmt6"


# Use HMMER Pfam database search to extract additional ORFs
hmmscan --cpu 24 --domtblout $shortName"pfam.domtblout" ${db_inpath}/Pfam-A.hmm ${transdecoder_dir}/$shortName".transdecoder.pep"

# Predict likely coding regions and store in output file
TransDecoder.Predict -t ${data_inpath}/$shortName".fa" --retain_pfam_hits $shortName".pfam.domtblout" \
                     --retain_blastp_hits $shortName".blastp.outfmt6"



end=`date +%s`
runtime=$((end-start))
echo $runtime
