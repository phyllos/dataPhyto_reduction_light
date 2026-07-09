#!/bin/bash

#SBATCH --job-name=tr_sfx_light_job01
#SBATCH --time=12:00:00
#SBATCH --mem=120G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --partition=bioxfel
#SBATCH --output=/home/uwm/gengh/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/logs/%x-%j.out
#SBATCH --error=/home/uwm/gengh/Data-cxfel/phytochrome2026/dataPhyto_reduction_light/logs/%x-%j.err

set -euo pipefail

echo "========================================"
echo "Job started at $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Running on $(hostname)"
echo "========================================"

source /home/uwm/gengh/Data-cxfel/conda/etc/profile.d/conda.sh
conda activate phytochrome_conda_env

PROJECT_DIR=/home/uwm/gengh/Data-cxfel/phytochrome2026/dataPhyto_reduction_light
cd "$PROJECT_DIR"



echo "Working directory: $(pwd)"
echo "Conda env: ${CONDA_DEFAULT_ENV:-NA}"
echo "Python: $(which python)"
echo "========================================"
echo

## Full C1 -> C13 run:
python ./driver_light_01.py --config ./inputs/config_light_upto1ps_01.yaml

## For a faster first test from existing C7 output, comment the line above and use:
## python ./driver_light.py --config ./inputs/config_light_from_c7.yaml

echo
echo "========================================"
echo "Job finished at $(date)"
echo "========================================"
