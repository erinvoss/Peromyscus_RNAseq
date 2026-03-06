#!/bin/bash

#SBATCH --job-name=trinity_peca
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --partition=savio3_bigmem
#SBATCH --account=fc_nachman
#SBATCH --output=slurmout/trinity_peca_job_%A.out
#SBATCH --error=slurmout/trinity_peca_job_%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erinvoss@berkeley.edu

cd /global/scratch/users/erinvoss/Peromyscus-RNA

start=`date +%s`
echo $HOSTNAME

# Re-run this script, substituting appropriate samples.txt file to assemble additional transcriptomes

inpath="/global/scratch/users/erinvoss/Peromyscus-RNA/01-FastP_Preproc" 
outpath="/global/scratch/users/erinvoss/Peromyscus-RNA/02-Trinity-Assembly/PECA-Trinity-Output"
sample_file="/global/scratch/users/erinvoss/Peromyscus-RNA/02-Trinity-Assembly/PECA_transcriptome_samples.txt"

# Initialize arrays
LEFT_FILES=()
RIGHT_FILES=()

# Read in cleaned fastq filenames 
while read -r SAMPLE; do
        [[ -z "$SAMPLE" ]] && continue
        [[ "$SAMPLE" =~ ^# ]] && continue
        
        R1="${inpath}/${SAMPLE}/${SAMPLE}_R1_cleaned.fastq.gz"
        R2="${inpath}/${SAMPLE}/${SAMPLE}_R2_cleaned.fastq.gz"

         # Safety check to make sure both files exist 
         for f in "$R1" "$R2"; do
                [[ -f "$f" ]] || { echo "Missing $f"; exit 1; }
        done

    # Define right and left files
    LEFT_FILES+=("$R1")
    RIGHT_FILES+=("$R2")
done < "$sample_file"

# Convert arrays to comma-separated lists 
LEFT=$(IFS=,; echo "${LEFT_FILES[*]}")
RIGHT=$(IFS=,; echo "${RIGHT_FILES[*]}")

# Run Trinity to assemble transcriptome
Trinity \
  --seqType fq \
  --max_memory 364G \
  --left  "$LEFT" \
  --right "$RIGHT" \
  --CPU 32 \
  --output "$outpath"

end=`date +%s`
runtime=$((end-start))
echo $runtime
