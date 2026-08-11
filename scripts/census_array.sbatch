#!/bin/bash
#SBATCH --job-name=psor_census
#SBATCH --partition=dept_cpu
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --output=logs/census_%A_%a.out
#SBATCH --error=logs/census_%A_%a.err
#SBATCH --mail-user=sba50@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
# NOTE: without ARRAY_TASKS in --mail-type, SLURM mails for the array AS A WHOLE
# (3 messages total), not once per task. Add ARRAY_TASKS only if you want per-task mail.

# recount3 blood-psoriasis census, sharded 4 ways.
# Deliberately modest: this is a shared cluster and the work is network-bound on
# an external server (recount3/duffel.rail.bio), so more shards mostly means more
# concurrent load on someone else's host, not proportionally faster completion.
# SLURM owns these processes, so they survive the submitting session ending --
# background procs spawned from a transient ssh job shell do not.
RS=/net/dali/home/mscbio/sba50/miniconda3/envs/scissor-r/bin/Rscript
cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p census_par logs

# Task 0 builds the shared project list; the rest wait for it (script handles this).
$RS code/census_blood_par.R "$SLURM_ARRAY_TASK_ID" 4 census_par
