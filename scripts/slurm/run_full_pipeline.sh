#!/bin/bash
#SBATCH --job-name=xenium_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=96:00:00
#SBATCH --output=logs/full_pipeline_%j.out
#SBATCH --error=logs/full_pipeline_%j.err
#SBATCH --qos=your_qos
#SBATCH --account=your_account

# Complete Xenium Analysis Pipeline for HiPerGator
# This script runs all analysis steps sequentially

echo "=========================================="
echo "Xenium Complete Analysis Pipeline"
echo "=========================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="

# Load modules
module purge
module load conda

# Activate environment
source activate xenium_analysis

# Create necessary directories
mkdir -p data/raw data/processed figures logs notebooks/executed

# Set environment variables
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Function to run notebook and check status
run_notebook() {
    local notebook=$1
    local name=$2
    
    echo ""
    echo "=========================================="
    echo "Running: $name"
    echo "Notebook: $notebook"
    echo "Started at: $(date)"
    echo "=========================================="
    
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=86400 \
        --output-dir=notebooks/executed \
        "notebooks/$notebook"
    
    if [ $? -eq 0 ]; then
        echo "SUCCESS: $name completed at $(date)"
    else
        echo "ERROR: $name failed with exit code $?"
        exit 1
    fi
}

# Run pipeline steps
run_notebook "01_preprocessing.ipynb" "Preprocessing"
run_notebook "02_phenotyping.ipynb" "Phenotyping"
run_notebook "03_spatial_analysis.ipynb" "Spatial Analysis"
run_notebook "04_group_comparisons.ipynb" "Group Comparisons"
run_notebook "05_tissue_comparisons.ipynb" "Tissue Comparisons"
run_notebook "06_xenium_phenocycler_integration.ipynb" "Phenocycler Integration"

echo ""
echo "=========================================="
echo "All analysis steps completed successfully!"
echo "Job finished at: $(date)"
echo "=========================================="
