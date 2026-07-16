# Compute Context

## Target
- Machine / cluster: **CSB cluster** (`cluster.csb.pitt.edu`); compute node allocated per session (this session: **n003**, 64 cores / 93 GB RAM).
- Access model: SSH to head node → `srun --pty -p dept_cpu -w <node> /bin/bash` to grab a node → register that node as a Claude Science SSH target (`ssh:<node>`). Node name changes each session; re-register each time. GlobalProtect VPN required off-campus.
- Module system: none used — conda sourced directly from `/net/dali/home/mscbio/sba50/miniconda3`.
- Submission method: `sbatch` via `call_command` on the registered node (the reliable route; `submit_job` auto-harvest is blocked by a probe timeout on CSB).
- Partitions: `dept_cpu` (used for the Scissor full-census run), `dept_gpu` (+`--gres=gpu:1`), `big_memory`, `any_cpu`/`any_gpu` (preemptable).
- GPU note: the Scissor solve is coordinate-descent elastic net (glmnet, CPU-only, no GPU backend) — CPU partition is correct. GPU is only relevant to the future sequence-level (Sei-LLRA) arm.
- Working folder: `~/Soahum/Project/psoriasis/psoriasis-1-bulk`. Heavy I/O scratch: `/scr/${SLURM_JOB_ID}`.
- Notification email: sba50@pitt.edu (`--mail-user=sba50@pitt.edu --mail-type=BEGIN,END,FAIL`).

## Environment (reproducibility-critical)
- Env name: **scissor-r** (`conda env create -f environment.yml && conda activate scissor-r`).
- Interpreter on node: `/net/dali/home/mscbio/sba50/miniconda3/envs/scissor-r/bin/Rscript`.
- Pins (from `environment.yml`):
```
r-base=4.5
r-seurat                        # Seurat 5
r-matrix
r-glmnet                        # network-regularized elastic net solver
r-rspectra
r-ggplot2 / r-patchwork / r-ggrepel
r-r.utils / r-data.table
bioconductor-preprocesscore     # quantile normalization
bioconductor-org.hs.eg.db       # Ensembl -> symbol
bioconductor-edger / bioconductor-limma
bioconductor-summarizedexperiment
bioconductor-recount3           # bulk anchor retrieval
bioconductor-geoquery           # GSE173706 download
```
- Local (laptop) bulk meta-analysis env: `psoriasis-r` (R 4.5.3; recount3/DESeq2/edgeR/limma) — used for the k=4 PP-vs-NN work.
- Solver choice: compiled Scissor **APML1 is NOT used**; replaced by a pure-R glmnet reimplementation of the network-regularized elastic net (`code/scissor_glmnet_solver.R`). Validated (synthetic recovery cor=0.72); canonical-solver cross-check still staged.
- Seeds: set inside the run scripts; Scissor alpha=0.40 selected by the tuning grid.
