#!/bin/bash
# Local (non-SLURM) pipeline driver for Xenium analysis.
# Runs preprocessing (01), phenotyping (02), and spatial (03) per sample.
# Integration (05) is NOT part of the default flow — per-sample analysis is
# the intended workflow. Designed to run inside the currently allocated
# HiPerGator VSCode tunnel session — DO NOT sbatch this script.
#
# Usage (from a terminal attached to the tunnel):
#   bash scripts/run_local_pipeline.sh              # default stages, all samples
#   bash scripts/run_local_pipeline.sh 01 0041323   # one stage, one sample
#   STAGES="01 02" SAMPLES="0041323" \
#     bash scripts/run_local_pipeline.sh             # env override
#   STAGES="05" bash scripts/run_local_pipeline.sh   # opt-in to integration
#
# Logs:  logs/local_<stage>_<sample>_<timestamp>.{out,err}
# Executed notebooks:  notebooks/executed/<sample>/<stage>_<sample>.ipynb
# Integrated notebook: notebooks/executed/integrated_<tag>/05_tissue_comparisons_<tag>.ipynb

set -euo pipefail

PROJECT_ROOT="/blue/maigan/smith6jt/Xenium_Analysis"
cd "$PROJECT_ROOT"

source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis

NCPU=$(nproc)
export OMP_NUM_THREADS="$NCPU"
export MKL_NUM_THREADS="$NCPU"
export NUMEXPR_MAX_THREADS="$NCPU"
export OPENBLAS_NUM_THREADS="$NCPU"
# scVI / PyTorch: use the single B200 GPU by default.
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"

# CUDA toolkit env for rapids_singlecell / cuml / cupy / cudf.
# See scripts/env_rapids.sh for the full rationale (CLAUDE.md §7).
source "${PROJECT_ROOT}/scripts/env_rapids.sh"

SAMPLES_DEFAULT="0041323 0041326"
SAMPLES="${SAMPLES:-$SAMPLES_DEFAULT}"

# Allow positional single-stage / single-sample invocation.
if [ $# -ge 1 ]; then
    STAGES="$1"
fi
if [ $# -ge 2 ]; then
    SAMPLES="$2"
fi
STAGES="${STAGES:-01 02 03}"

mkdir -p logs

run_per_sample () {
    local stage="$1"
    local sample="$2"
    local nb
    case "$stage" in
        01) nb="notebooks/01_preprocessing_v2.ipynb" ;;
        02) nb="notebooks/02_phenotyping.ipynb" ;;
        03) nb="notebooks/03_spatial_analysis.ipynb" ;;
        *)  echo "unknown per-sample stage: $stage" >&2; return 2 ;;
    esac

    local ts
    ts=$(date +%Y%m%d_%H%M%S)
    local outdir="notebooks/executed/${sample}"
    local figdir="figures/${sample}"
    mkdir -p "$outdir" "${figdir}/01_preprocessing" "${figdir}/02_phenotyping" "${figdir}/03_spatial_analysis"

    local out="logs/local_${stage}_${sample}_${ts}.out"
    local err="logs/local_${stage}_${sample}_${ts}.err"

    echo "[$(date)] START stage=$stage sample=$sample"
    echo "  notebook: $nb"
    echo "  logs:     $out"
    echo "            $err"

    XENIUM_SAMPLE_ID="$sample" \
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=86400 \
        --output-dir="$outdir" \
        --output="${stage}_${sample}.ipynb" \
        "$nb" >"$out" 2>"$err"

    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "[$(date)] FAIL stage=$stage sample=$sample (rc=$rc)"
        echo "  tail of $err:"
        tail -30 "$err" || true
        return "$rc"
    fi
    echo "[$(date)] OK   stage=$stage sample=$sample"
}

run_integration () {
    local sample_list="$1"   # comma-separated
    local tag
    tag=$(echo "$sample_list" | tr ',' '_')

    local ts
    ts=$(date +%Y%m%d_%H%M%S)
    local outdir="notebooks/executed/integrated_${tag}"
    local figdir="figures/integrated_${tag}"
    mkdir -p "$outdir" "$figdir" "data/processed/integrated"

    local out="logs/local_05_${tag}_${ts}.out"
    local err="logs/local_05_${tag}_${ts}.err"

    echo "[$(date)] START stage=05 samples=$sample_list"
    echo "  logs: $out"
    echo "        $err"

    XENIUM_SAMPLE_IDS="$sample_list" \
    XENIUM_FIGURES_DIR="$figdir" \
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=86400 \
        --output-dir="$outdir" \
        --output="05_tissue_comparisons_${tag}.ipynb" \
        notebooks/05_tissue_comparisons.ipynb >"$out" 2>"$err"

    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "[$(date)] FAIL stage=05 (rc=$rc)"
        tail -30 "$err" || true
        return "$rc"
    fi
    echo "[$(date)] OK   stage=05"
}

for stage in $STAGES; do
    case "$stage" in
        01|02|03)
            for sample in $SAMPLES; do
                run_per_sample "$stage" "$sample"
            done
            ;;
        05)
            csv=$(echo "$SAMPLES" | tr ' ' ',')
            run_integration "$csv"
            ;;
        *)
            echo "Unknown stage: $stage (expected 02, 03, or 05)" >&2
            exit 2
            ;;
    esac
done

echo "[$(date)] Pipeline finished. stages=[$STAGES] samples=[$SAMPLES]"
