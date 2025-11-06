#!/bin/bash
#SBATCH --job-name=xenium_spatial
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=48:00:00
#SBATCH --output=logs/spatial_%j.out
#SBATCH --error=logs/spatial_%j.err
#SBATCH --qos=your_qos
#SBATCH --account=your_account

# Xenium Spatial Analysis Pipeline for HiPerGator
# This script runs comprehensive spatial analysis with Squidpy

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load modules
module purge
module load conda

# Activate environment
source activate xenium_analysis

# Set environment variables for parallel processing
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run spatial analysis notebook
echo "Starting spatial analysis..."
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=172800 \
    --output-dir=notebooks/executed \
    notebooks/03_spatial_analysis.ipynb

if [ $? -eq 0 ]; then
    echo "Spatial analysis completed successfully"
else
    echo "Spatial analysis failed with exit code $?"
    exit 1
fi

echo "Job finished at: $(date)"
