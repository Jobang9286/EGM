#!/bin/bash
#SBATCH --job-name=EGM_preprocess
#SBATCH --cpus-per-task=48
#SBATCH --mem=1000G
#SBATCH --time=06:00:00
#SBATCH --output=EGM_data/preprocess_%j.out
#SBATCH --error=EGM_data/preprocess_%j.err
echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate senepy

# Install scrublet if not present
pip install scrublet --quiet 2>/dev/null

# Run preprocessing
cd EGM_data
python preprocess_EGM.py

echo "Job finished: $(date)"
