#!/bin/bash
# Jupyter kernel launcher for the xenium_analysis env with the rapids env fix
# pre-applied. Referenced by ~/.local/share/jupyter/kernels/xenium-rapids/kernel.json
# so that LD_PRELOAD is set BEFORE the python process starts (it cannot be set
# from inside an already-running interpreter).
#
# Args: ipykernel_launcher passes through the connection-file flags that Jupyter
# substitutes for {connection_file}.

set -euo pipefail

source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis

source /blue/maigan/smith6jt/Xenium_Analysis/scripts/env_rapids.sh

exec python -m ipykernel_launcher "$@"
