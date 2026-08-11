#!/bin/bash
#SBATCH --job-name=scissor-full-census
#SBATCH --partition=big_memory
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --output=logs/full_census_%j.out
#SBATCH --error=logs/full_census_%j.err
#
# Full-census Scissor-on-gradient run (all 89,058 QC-passed cells).
# See code/run_full_census_cluster.R for the pipeline; it is checkpointed
# per-stage under results_full/, so if this job is killed (wall-time,
# `scontrol requeue`, node failure) simply `sbatch` it again — it resumes
# at the first incomplete stage instead of restarting from scratch.
#
# Submit from the repo root:
#   sbatch code/run_full_census.sbatch

set -euo pipefail
mkdir -p logs results_full

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scissor-r

Rscript code/run_full_census_cluster.R
