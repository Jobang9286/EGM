#!/bin/bash
#SBATCH --job-name=hdwgcna
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=EGM_data/hdwgcna_pipeline/hdwgcna_%j.out
#SBATCH --error=EGM_data/hdwgcna_pipeline/hdwgcna_%j.err
cd EGM_data/hdwgcna_pipeline
source $(conda info --base)/etc/profile.d/conda.sh
conda activate senepy
python -u run_hdwgcna.py
