#!/bin/bash
#SBATCH --job-name=psor_census
#SBATCH --partition=dept_cpu
#SBATCH --array=0-15
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=logs/census_%A_%a.out
#SBATCH --error=logs/census_%A_%a.err

# recount3 blood-psoriasis census, sharded 16 ways.
# SLURM owns these processes, so they survive the submitting session ending --
# background procs spawned from a transient ssh job shell do not.
RS=/net/dali/home/mscbio/sba50/miniconda3/envs/scissor-r/bin/Rscript
cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p census_par logs

# Task 0 builds the shared project list; the rest wait for it (script handles this).
$RS code/census_blood_par.R "$SLURM_ARRAY_TASK_ID" 16 census_par
