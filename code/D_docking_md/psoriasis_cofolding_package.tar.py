# Auto-extracted generating script
# Produces: psoriasis_cofolding_package.tar.gz
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): validation.json, domain_sequences.json, docking_ranked.csv, vyi_refH.sdf, RORC_5APH.pdb, screening_library.csv
# Source artifact version: 13144510-cd08-41b7-88cd-5d55f66391bb
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import json, os, tarfile, io
from rdkit import Chem

os.makedirs("work/cluster_packages/02_cofolding/yaml", exist_ok=True)
os.makedirs("work/cluster_packages/02_cofolding/scripts", exist_ok=True)
os.makedirs("work/cluster_packages/02_cofolding/sequences", exist_ok=True)
os.makedirs("work/cluster_packages/02_cofolding/logs", exist_ok=True)
os.makedirs("work/cluster_packages/02_cofolding/out", exist_ok=True)

doms = json.load(open("domain_sequences.json"))

# clean RORC TEV tag
rorc = doms["RORC"]["domain_seq"]
if "ENLYFQG" in rorc[:15]:
    rorc = rorc[rorc.index("ENLYFQG")+7:]
    rorc = rorc.lstrip("AS")
doms["RORC"]["domain_seq_clean"] = rorc

# Save domain sequences as FASTA reference
with open("work/cluster_packages/02_cofolding/sequences/domain_constructs.fasta", "w") as f:
    f.write(f">STAT3_core_SH2_6NJS_chainA {doms['STAT3']['len']}aa\n{doms['STAT3']['domain_seq']}\n")
    f.write(f">RORC_LBD_5APH_clean {len(doms['RORC']['domain_seq_clean'])}aa\n{doms['RORC']['domain_seq_clean']}\n")
    f.write(f">JAK3_kinase_5LWM_chainA {doms['JAK3']['len']}aa\n{doms['JAK3']['domain_seq']}\n")

# template YAML for STAT3 and JAK3
for sym, key in [("STAT3", "domain_seq"), ("JAK3", "domain_seq")]:
    seq = doms[sym][key]
    tmpl = f"""# Boltz-2 co-folding + affinity template for {sym}
# Fill in top docking-hit SMILES (from package 01 results) — one YAML per ligand.
version: 1
sequences:
  - protein:
      id: A
      sequence: {seq}
  - ligand:
      id: L
      smiles: 'REPLACE_WITH_HIT_SMILES'
properties:
  - affinity:
      binder: L
"""
    open(f"work/cluster_packages/02_cofolding/yaml/_TEMPLATE_{sym}.yaml", "w").write(tmpl)

import pandas as pd
valid = pd.read_csv("docking_ranked.csv")
top = valid.head(10)[["chembl_id", "smiles", "affinity", "mw", "qed"]].copy()

vyi = Chem.MolFromMolFile("vyi_refH.sdf")
vyi_smi = Chem.MolToSmiles(vyi) if vyi else None

def n_heavy(smi):
    m = Chem.MolFromSmiles(smi)
    return m.GetNumHeavyAtoms() if m else 999

hits = []
if vyi_smi:
    hits.append(("VYI_crystal_posctrl", vyi_smi, "positive_control"))
for _, r in top.iterrows():
    hits.append((r.chembl_id, r.smiles, f"dock{r.affinity:.1f}"))

written = 0
for name, smi, tag in hits:
    nh = n_heavy(smi)
    if nh > 128:
        print(f"  skip {name}: {nh} heavy atoms > 128 (Boltz affinity cap)")
        continue
    yaml = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {rorc}
  - ligand:
      id: L
      smiles: '{smi}'
properties:
  - affinity:
      binder: L
"""
    open(f"work/cluster_packages/02_cofolding/yaml/RORC_{name}.yaml", "w").write(yaml)
    written += 1
print(f"\nwrote {written} Boltz YAML inputs for RORC co-folding")

# run_boltz.slurm
slurm_content = """#!/bin/bash
#SBATCH --job-name=psor_boltz
#SBATCH --array=0-10
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/boltz_%A_%a.out
#SBATCH --error=logs/boltz_%A_%a.err
# ---------------------------------------------------------------------------
# Boltz-2 co-folding + affinity cross-check on top docking hits (GPU).
# One array task per YAML input. Each predicts the protein-ligand complex
# structure AND a binding-affinity estimate (affinity head).
#
# EDIT for your cluster: GPU partition/gres string, conda/module activation.
# The MSA step queries api.colabfold.com (needs outbound network on the node,
# or pre-compute .a3m files and reference them under msa: in each YAML).
# ---------------------------------------------------------------------------
set -euo pipefail
mkdir -p logs out

# --- Environment (EDIT) ---
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null && conda activate boltz-env || true
# pip install boltz   # if not yet installed

# pick the YAML for this array index
YAMLS=(yaml/RORC_*.yaml)
YAML=${YAMLS[$SLURM_ARRAY_TASK_ID]}
echo "Predicting: $YAML"

boltz predict "$YAML" \\
    --use_msa_server \\
    --out_dir out/ \\
    --recycling_steps 3 \\
    --diffusion_samples 5 \\
    --output_format pdb
    # add --no_kernels if cuequivariance_ops_torch ImportError

# Outputs land in out/boltz_results_<name>/predictions/<name>/
#   confidence_<name>_model_0.json   (iptm, complex_plddt, confidence_score)
#   affinity_<name>.json             (affinity_pred_value, affinity_probability_binary)
#   <name>_model_0.pdb               (predicted complex)
"""
open("work/cluster_packages/02_cofolding/scripts/run_boltz.slurm", "w").write(slurm_content)

# aggregate_cofolding.py
aggregate_content = """#!/usr/bin/env python
\"\"\"
Aggregate Boltz-2 co-folding + affinity outputs into one ranked table and
cross-check against the docking scores.

Usage:
    python aggregate_cofolding.py --out_dir out --docking docking_ranked.csv \\
        --result cofolding_crosscheck.csv
\"\"\"
import os, glob, json, argparse, csv

def find(pred_dir, pat):
    hits=glob.glob(os.path.join(pred_dir,pat))
    return hits[0] if hits else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out_dir",default="out")
    ap.add_argument("--docking",default=None,help="optional docking_ranked.csv to merge")
    ap.add_argument("--result",default="cofolding_crosscheck.csv")
    args=ap.parse_args()

    dock={}
    if args.docking and os.path.exists(args.docking):
        for r in csv.DictReader(open(args.docking)):
            dock[r["ligand"]]=r

    rows=[]
    for rd in sorted(glob.glob(f"{args.out_dir}/boltz_results_*")):
        name=os.path.basename(rd).replace("boltz_results_","")
        pred=os.path.join(rd,"predictions",name)
        conf=find(pred,"confidence_*_model_0.json")
        aff =find(pred,"affinity_*.json")
        row={"ligand":name.replace("RORC_","")}
        if conf:
            c=json.load(open(conf))
            row["iptm"]=c.get("iptm"); row["complex_plddt"]=c.get("complex_plddt")
            row["confidence_score"]=c.get("confidence_score")
        if aff:
            a=json.load(open(aff))
            row["affinity_pred_log10_uM"]=a.get("affinity_pred_value")
            row["affinity_prob_binder"]=a.get("affinity_probability_binary")
        if row["ligand"] in dock:
            row["vina_affinity"]=dock[row["ligand"]].get("affinity")
        rows.append(row)

    # rank by Boltz binder probability
    rows.sort(key=lambda r: -(r.get("affinity_prob_binder") or 0))
    if rows:
        cols=["ligand","affinity_prob_binder","affinity_pred_log10_uM","iptm",
              "complex_plddt","confidence_score","vina_affinity"]
        with open(args.result,"w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=cols,extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        print(f"{len(rows)} complexes -> {args.result}")
        for r in rows[:5]:
            print(f"  {r['ligand']}: P(binder)={r.get('affinity_prob_binder')} "
                  f"ipTM={r.get('iptm')} vina={r.get('vina_affinity')}")
    else:
        print("no Boltz outputs found under", args.out_dir)

if __name__=="__main__":
    main()
"""
open("work/cluster_packages/02_cofolding/scripts/aggregate_cofolding.py", "w").write(aggregate_content)

# README.md
readme_content = """# Co-folding + affinity cross-check — GPU cluster package

Orthogonal, structure-based validation of the top docking hits using **Boltz-2**
(open-weights AF3-class co-folder). For each protein–ligand pair Boltz-2 predicts
the **complex structure** (with ipTM/pLDDT confidence) and a **binding-affinity
estimate** via its affinity head — an independent check that does not share Vina's
scoring assumptions.

## Contents