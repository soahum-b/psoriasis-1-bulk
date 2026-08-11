# scripts/

Every runnable shell script in the project, in one place. These are **copies** of
the originals under `code/` — kept in both locations so nothing that references a
`code/…` path breaks (the project follows an additive-only rule). If you edit one,
sync the other.

| script | submits / does | job name | notes |
|---|---|---|---|
| `census_array.sbatch` | recount3 blood-psoriasis census, `code/census_blood_par.R` sharded 4 ways | `psor_census` | 4-task array, `dept_cpu`, one node, 12 h walltime. Emails `sba50@pitt.edu` on BEGIN/END/FAIL (3 mails for the array as a whole, not per task). Completed run: array **56919691**. |
| `run_full_census.sbatch` | full-census Scissor run, `code/run_full_census_cluster.R` | `scissor-full-census` | Completed run: job **56882314** → `results_full/`. CPU-only (glmnet coordinate descent has no GPU path). |
| `regen_fig4.sh` | regenerates Figure 4 from the full-census object | `psor_fig4` | Small single job. |
| `pipeline_scripts/setup.sh` | one-command conda environment setup for the bulk meta-analysis | — | Not a SLURM script; run interactively. |

## Cluster notes (CSB)

Non-interactive shells lack some tooling on PATH:

    export PATH=/opt/slurm/bin:$PATH      # sbatch, squeue, sacct
    bash -lc 'module avail'               # module is a login-shell function

R lives in conda envs, not a module: `~/miniconda3/envs/scissor-r` (R 4.5.3,
Seurat + recount3) and `~/miniconda3/envs/r441_env` (R 4.4.1, DESeq2 + recount3).
There is no `psoriasis-r` env on the cluster — that name is local-sandbox only.

**Submit long work with `sbatch`, not background processes.** Detaching with
`setsid nohup … &` from a transient ssh session is unreliable: the children get
reaped when that session ends (observed — 16 spawned workers left 1 alive).
