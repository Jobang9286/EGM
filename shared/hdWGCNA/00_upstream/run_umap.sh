#!/bin/bash
#SBATCH --job-name=EGM_UMAP
#SBATCH --cpus-per-task=48
#SBATCH --mem=1000G
#SBATCH --time=04:00:00
#SBATCH --output=EGM_data/senepy_pipeline/umap_%j.out
#SBATCH --error=EGM_data/senepy_pipeline/umap_%j.err
echo "Job started: $(date)"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate senepy
pip install harmonypy --quiet 2>/dev/null
cd EGM_data/senepy_pipeline
python -u make_umap.py
echo "Job finished: $(date)"
