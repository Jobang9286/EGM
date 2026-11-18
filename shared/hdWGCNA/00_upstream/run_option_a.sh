#!/bin/bash
#SBATCH --job-name=EGM_OptA
#SBATCH --cpus-per-task=48
#SBATCH --mem=512G
#SBATCH --time=02:00:00
#SBATCH --output=EGM_data/senepy_pipeline/optA_%j.out
#SBATCH --error=EGM_data/senepy_pipeline/optA_%j.err
echo "Job started: $(date)"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate senepy
python3 -u EGM_data/senepy_pipeline/option_a_stratify.py
echo "Job finished: $(date)"
