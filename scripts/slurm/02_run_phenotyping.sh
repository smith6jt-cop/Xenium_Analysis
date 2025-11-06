#!/bin/bash
#SBATCH --job-name=xenium_phenotype
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/phenotype_%j.out
#SBATCH --error=logs/phenotype_%j.err
#SBATCH --qos=your_qos
#SBATCH --account=your_account

# Xenium Phenotyping Pipeline for HiPerGator
# This script runs cell type annotation with scVI

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load modules
module purge
module load conda

# Activate environment
source activate xenium_analysis

# Set environment variables
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run phenotyping notebook
echo "Starting phenotyping analysis..."
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=86400 \
    --output-dir=notebooks/executed \
    notebooks/02_phenotyping.ipynb

if [ $? -eq 0 ]; then
    echo "Phenotyping completed successfully"
else
    echo "Phenotyping failed with exit code $?"
    exit 1
fi

echo "Job finished at: $(date)"
