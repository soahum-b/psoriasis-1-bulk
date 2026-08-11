#!/bin/bash
#SBATCH -J psor_fig4
#SBATCH -p dept_cpu
#SBATCH -t 00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH -o psor_fig4_%j.out
eval "$(conda shell.bash hook)"
conda activate scissor-r
cd /net/dali/home/mscbio/sba50/Soahum/Project/psoriasis/psoriasis-1-bulk
Rscript code/regen_fig4_fullcensus.R
echo "EXIT $?"
