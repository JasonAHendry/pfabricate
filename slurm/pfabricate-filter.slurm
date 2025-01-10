#!/bin/bash -l
#SBATCH --job-name=nbar
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

# Fixed width for chomrosome
CHROM=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

pfabricate filter \
  -i {in_vcf_template} \
  -o {out_vcf} \
  -s {sample_str} \


