#!/bin/bash

#SBATCH --job-name=phytochrome_light_job03
#SBATCH --time=12:00:00
#SBATCH --mem=400G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=cxfel
#SBATCH --nodelist=execute-4005
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

echo "========================================"
echo "Job started at $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Running on $(hostname)"
echo "========================================"

PROJECT_DIR="${PROJECT_DIR:-$HOME/Data-cxfel/phytochrome2026/dataPhyto_reduction_light}"
CONDA_BASE="${CONDA_BASE:-$HOME/Data-cxfel/conda}"
CONDA_ENV="${CONDA_ENV:-phytochrome_conda_env}"
CONFIG_FILE="${CONFIG_FILE:-inputs/config_light_upto1ps_02.yaml}"

cd "$PROJECT_DIR"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"




echo "Working directory: $(pwd)"
echo "Conda env: ${CONDA_DEFAULT_ENV:-NA}"
echo "Python: $(which python)"
echo "Config file: $CONFIG_FILE"
echo "========================================"
echo

## Full C1 -> C13 run:
python ./driver_light_01.py --config "$CONFIG_FILE"

## For a faster first test from existing C7 output, comment the line above and use:
## python ./driver_light.py --config ./inputs/config_light_from_c7.yaml

echo
echo "========================================"
echo "Job finished at $(date)"
echo "========================================"
