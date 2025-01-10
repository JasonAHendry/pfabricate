#!/bin/bash -l
#SBATCH --job-name=pfabfilt
#SBATCH --output=logs/filter-jid%A-%a.out
#SBATCH --error=logs/filter-jid%A-%a.err
#SBATCH --chdir=./
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array {array_str}
#SBATCH --mem=16GB
#SBATCH --time=23:00:00

# LOAD ANACONDA
module load anaconda/3/.2023.03

# ACTIVATE ENVIRONMENT
conda activate pfabricate

# Prepare variables
CHROM=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
INPUT_VCF={in_vcf_template}
OUTPUT_VCF={out_vcf}

echo "Processing chromosome: " $CHROM
echo "Input VCF: " $INPUT_VCF
echo "Output VCF: " $OUTPUT_VCF

pfabricate filter \
  -i $INPUT_VCF \
  -o $OUTPUT_VCF \
  -s {sample_str} \


