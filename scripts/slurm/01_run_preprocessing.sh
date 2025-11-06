#!/bin/bash
#SBATCH --job-name=xenium_preprocess
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/preprocess_%j.out
#SBATCH --error=logs/preprocess_%j.err
#SBATCH --qos=your_qos
#SBATCH --account=your_account

# Xenium Preprocessing Pipeline for HiPerGator
# This script runs the preprocessing notebook on HiPerGator cluster

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load required modules
module purge
module load conda

# Activate conda environment
source activate xenium_analysis

# Create output directories
mkdir -p data/raw data/processed figures/01_preprocessing logs

# Set environment variables
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run preprocessing notebook
echo "Starting preprocessing..."
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=86400 \
    --output-dir=notebooks/executed \
    notebooks/01_preprocessing.ipynb

# Check exit status
if [ $? -eq 0 ]; then
    echo "Preprocessing completed successfully"
else
    echo "Preprocessing failed with exit code $?"
    exit 1
fi

echo "Job finished at: $(date)"
