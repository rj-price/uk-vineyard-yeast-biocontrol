#!/bin/bash
#SBATCH --job-name=allelic
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=analysis/allelic/logs/allelic_%j.log
set -euo pipefail
cd /mnt/shared/scratch/jnprice/private/biocontrol_yeast
conda run -n data_science python analysis/scripts/15_allelic.py
