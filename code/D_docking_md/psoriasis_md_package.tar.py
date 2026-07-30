# Auto-extracted generating script
# Produces: psoriasis_md_package.tar.gz
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): RORC_apo.pdb, equilibrated.pdb, system.xml, state.xml
# Source artifact version: 08f0a982-879f-409a-8caa-c242b9798fac
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import os, shutil, tarfile

# Create directory structure
os.makedirs("work/cluster_packages/03_md/systems/RORC_crystal_VYI", exist_ok=True)
os.makedirs("work/cluster_packages/03_md/systems/RORC_tophit_CHEMBL3314024", exist_ok=True)
os.makedirs("work/cluster_packages/03_md/scripts", exist_ok=True)
os.makedirs("work/cluster_packages/03_md/logs", exist_ok=True)

# Copy equilibrated system files for RORC_tophit_CHEMBL3314024
shutil.copy("equilibrated.pdb", "work/cluster_packages/03_md/systems/RORC_tophit_CHEMBL3314024/equilibrated.pdb")
shutil.copy("system.xml", "work/cluster_packages/03_md/systems/RORC_tophit_CHEMBL3314024/system.xml")
shutil.copy("state.xml", "work/cluster_packages/03_md/systems/RORC_tophit_CHEMBL3314024/state.xml")

# Write build_system.py
build_system_py = '''#!/usr/bin/env python
"""
Build a solvated, parameterized, energy-minimized MD system for a
protein-ligand complex and run NVT+NPT equilibration.

Amber ff14SB (protein) + OpenFF-2.x (ligand, via openmmforcefields) + TIP3P water.

Usage:
    python build_system.py --receptor RORC_apo.pdb --ligand_sdf ligand.sdf \\
        --out systems/RORC_tophit --pad 1.0 --ion_conc 0.15 \\
        --nvt_ps 100 --npt_ps 100

Outputs (in --out):
    system.xml       serialized OpenMM System (parameterized, solvated)
    equilibrated.pdb final equilibrated coordinates (production start state)
    state.xml        serialized State (positions + velocities + box)
    build_report.json
"""
import os, json, argparse, time

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--receptor",required=True,help="protein PDB (apo, no ligand/waters)")
    ap.add_argument("--ligand_sdf",required=True,help="ligand 3D SDF, correctly protonated, in the pocket")
    ap.add_argument("--out",required=True)
    ap.add_argument("--pad",type=float,default=1.0,help="solvent padding (nm)")
    ap.add_argument("--ion_conc",type=float,default=0.15,help="NaCl conc (M)")
    ap.add_argument("--nvt_ps",type=float,default=100.0)
    ap.add_argument("--npt_ps",type=float,default=100.0)
    ap.add_argument("--min_only",action="store_true",help="minimize only, skip equilibration (fast test)")
    ap.add_argument("--platform",default=None,help="CUDA|OpenCL|CPU (auto if unset)")
    args=ap.parse_args()
    os.makedirs(args.out,exist_ok=True)

    from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, Platform, XmlSerializer
    from openmm.app import PDBFile, Modeller, ForceField, PME, HBonds
    import pdbfixer
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SystemGenerator
    t0=time.time(); rep={"args":vars(args)}

    # 1) Fix receptor (missing atoms/residues, protonation at pH 7.4)
    fixer=pdbfixer.PDBFixer(filename=args.receptor)
    fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    prot=Modeller(fixer.topology, fixer.positions)
    rep["protein_atoms"]=prot.topology.getNumAtoms()

    # 2) Ligand -> OpenFF Molecule (parameterized on the fly by SystemGenerator)
    lig=Molecule.from_file(args.ligand_sdf)
    lig_top=lig.to_topology().to_openmm()
    lig_pos=lig.conformers[0].to_openmm()

    # 3) Combine protein + ligand
    prot.add(lig_top, lig_pos)

    # 4) SystemGenerator: ff14SB + OpenFF for the small molecule + TIP3P
    ff_kwargs=dict(constraints=HBonds, rigidWater=True, removeCMMotion=True,
                   hydrogenMass=1.5*unit.amu)
    sysgen=SystemGenerator(
        forcefields=["amber/ff14SB.xml","amber/tip3p_standard.xml"],
        small_molecule_forcefield="openff-2.1.0",
        molecules=[lig], forcefield_kwargs=ff_kwargs, cache="ligand_ff.json")

    # 5) Solvate + neutralize + salt
    modeller=Modeller(prot.topology, prot.positions)
    modeller.addSolvent(sysgen.forcefield, model="tip3p",
                        padding=args.pad*unit.nanometer,
                        ionicStrength=args.ion_conc*unit.molar,
                        neutralize=True)
    rep["total_atoms"]=modeller.topology.getNumAtoms()

    system=sysgen.create_system(modeller.topology)

    # 6) Integrator + platform
    integ=LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
    if args.platform:
        plat=Platform.getPlatformByName(args.platform); sim=app.Simulation(modeller.topology,system,integ,plat)
    else:
        sim=app.Simulation(modeller.topology,system,integ)
    sim.context.setPositions(modeller.positions)
    rep["platform"]=sim.context.getPlatform().getName()

    # 7) Minimize
    sim.minimizeEnergy(maxIterations=5000)
    rep["minimized"]=True

    if not args.min_only:
        # 8) NVT heat-up (restrained not shown for brevity; short NVT then NPT)
        sim.context.setVelocitiesToTemperature(300*unit.kelvin)
        sim.step(int(args.nvt_ps*1000/2))  # 2 fs steps
        rep["nvt_ps"]=args.nvt_ps
        # 9) NPT
        system.addForce(MonteCarloBarostat(1.0*unit.bar, 300*unit.kelvin, 25))
        sim.context.reinitialize(preserveState=True)
        sim.step(int(args.npt_ps*1000/2))
        rep["npt_ps"]=args.npt_ps

    # 10) Serialize equilibrated system + state
    with open(f"{args.out}/system.xml","w") as f: f.write(XmlSerializer.serialize(system))
    st=sim.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
    with open(f"{args.out}/state.xml","w") as f: f.write(XmlSerializer.serialize(st))
    with open(f"{args.out}/equilibrated.pdb","w") as f:
        PDBFile.writeFile(sim.topology, st.getPositions(), f)
    rep["wall_s"]=round(time.time()-t0,1)
    json.dump(rep, open(f"{args.out}/build_report.json","w"), indent=1)
    print(f"[{os.path.basename(args.out)}] built: {rep[\'total_atoms\']} atoms, "
          f"{rep[\'platform\']}, {rep[\'wall_s\']}s -> {args.out}")

if __name__=="__main__":
    main()
'''
open("work/cluster_packages/03_md/scripts/build_system.py", "w").write(build_system_py)

# Write run_production.py
run_production_py = '''#!/usr/bin/env python
"""
Run production MD from a pre-equilibrated system (system.xml + state.xml).

Usage:
    python run_production.py --sys systems/RORC_tophit_CHEMBL3314024 \\
        --ns 100 --report_ps 10 --platform CUDA

Reads system.xml + state.xml (positions+velocities+box) from --sys, runs NPT
production, writes a DCD trajectory, a state data log, and a final checkpoint.
"""
import os, argparse, time

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sys",required=True,help="dir containing system.xml + state.xml + equilibrated.pdb")
    ap.add_argument("--ns",type=float,default=100.0,help="production length (ns)")
    ap.add_argument("--dt_fs",type=float,default=2.0)
    ap.add_argument("--report_ps",type=float,default=10.0,help="trajectory/log interval (ps)")
    ap.add_argument("--temp",type=float,default=300.0)
    ap.add_argument("--platform",default="CUDA")
    ap.add_argument("--precision",default="mixed")
    args=ap.parse_args()

    from openmm import app, unit, LangevinMiddleIntegrator, Platform, XmlSerializer
    from openmm.app import PDBFile, Simulation, DCDReporter, StateDataReporter, CheckpointReporter

    toppdb=PDBFile(f"{args.sys}/equilibrated.pdb")
    with open(f"{args.sys}/system.xml") as f: system=XmlSerializer.deserialize(f.read())
    with open(f"{args.sys}/state.xml") as f: state=XmlSerializer.deserialize(f.read())

    integ=LangevinMiddleIntegrator(args.temp*unit.kelvin, 1.0/unit.picosecond, args.dt_fs*unit.femtoseconds)
    try:
        plat=Platform.getPlatformByName(args.platform)
        props={"Precision":args.precision} if args.platform=="CUDA" else {}
        sim=Simulation(toppdb.topology, system, integ, plat, props)
    except Exception:
        sim=Simulation(toppdb.topology, system, integ)  # auto platform fallback
    sim.context.setState(state)

    nsteps=int(args.ns*1e6/args.dt_fs)           # ns -> steps
    interval=int(args.report_ps*1000/args.dt_fs)  # ps -> steps
    out=args.sys
    sim.reporters.append(DCDReporter(f"{out}/production.dcd", interval))
    sim.reporters.append(StateDataReporter(f"{out}/production.log", interval,
        step=True, time=True, potentialEnergy=True, temperature=True, density=True,
        speed=True, remainingTime=True, totalSteps=nsteps))
    sim.reporters.append(CheckpointReporter(f"{out}/production.chk", interval*10))

    print(f"[{os.path.basename(out)}] production: {args.ns} ns = {nsteps} steps on "
          f"{sim.context.getPlatform().getName()}")
    t0=time.time(); sim.step(nsteps)
    with open(f"{out}/final_state.xml","w") as f: f.write(XmlSerializer.serialize(
        sim.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)))
    print(f"[{os.path.basename(out)}] DONE {args.ns} ns in {(time.time()-t0)/3600:.2f} h")

if __name__=="__main__":
    main()
'''
open("work/cluster_packages/03_md/scripts/run_production.py", "w").write(run_production_py)

# Write run_md.slurm
run_md_slurm = '''#!/bin/bash
#SBATCH --job-name=psor_md
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/md_%A_%a.out
#SBATCH --error=logs/md_%A_%a.err
#SBATCH --array=0-1
# ---------------------------------------------------------------------------
# Production MD from pre-equilibrated systems (one GPU per system).
# The systems in systems/ are already parameterized, solvated, minimized and
# equilibrated (NVT+NPT) — this script only runs production.
#
# EDIT: GPU gres string, conda/module activation, --ns length.
# ---------------------------------------------------------------------------
set -euo pipefail
mkdir -p logs

# --- Environment (EDIT) ---
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null && conda activate dock-md || true

SYSTEMS=(systems/RORC_crystal_VYI systems/RORC_tophit_CHEMBL3314024)
SYS=${SYSTEMS[$SLURM_ARRAY_TASK_ID]}

python scripts/run_production.py \\
    --sys "$SYS" \\
    --ns 100 \\
    --report_ps 10 \\
    --temp 300 \\
    --platform CUDA \\
    --precision mixed

# Analyze afterwards:
#   python scripts/analyze_traj.py --sys "$SYS"
'''
open("work/cluster_packages/03_md/scripts/run_md.slurm", "w").write(run_md_slurm)

# Write analyze_traj.py
analyze_traj_py = '''#!/usr/bin/env python
"""
Post-production trajectory analysis: protein backbone RMSD, per-residue RMSF,
ligand-RMSD (binding-pose stability), and protein-ligand contact count.

Usage:
    python analyze_traj.py --sys systems/RORC_tophit_CHEMBL3314024
Outputs (in --sys): analysis/{rmsd.csv,rmsf.csv,ligand_rmsd.csv,contacts.csv,md_analysis.png}
"""
import os, argparse
import numpy as np

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sys",required=True)
    ap.add_argument("--lig_resname",default="UNK",help="ligand residue name in topology")
    args=ap.parse_args()
    import mdtraj as md
    out=os.path.join(args.sys,"analysis"); os.makedirs(out,exist_ok=True)

    top=os.path.join(args.sys,"equilibrated.pdb")
    traj=md.load(os.path.join(args.sys,"production.dcd"), top=top)
    traj.image_molecules(inplace=True)
    prot=traj.atom_slice(traj.top.select("protein"))
    prot.superpose(prot, 0)

    # backbone RMSD
    bb=prot.top.select("backbone")
    rmsd=md.rmsd(prot, prot, 0, atom_indices=bb)*10  # nm->A
    np.savetxt(f"{out}/rmsd.csv", np.c_[traj.time/1000, rmsd], delimiter=",", header="time_ns,rmsd_A", comments="")

    # RMSF (per residue CA)
    ca=prot.top.select("name CA")
    rmsf=md.rmsf(prot, prot, 0, atom_indices=ca)*10
    resids=[prot.top.atom(i).residue.resSeq for i in ca]
    np.savetxt(f"{out}/rmsf.csv", np.c_[resids, rmsf], delimiter=",", header="resid,rmsf_A", comments="")

    # ligand RMSD (pose stability) — align on protein, measure ligand heavy-atom drift
    lig_sel=traj.top.select(f"resname {args.lig_resname} and not element H")
    figs=[("rmsd",traj.time/1000,rmsd,"Backbone RMSD (Å)")]
    if len(lig_sel)>0:
        traj.superpose(traj, 0, atom_indices=traj.top.select("protein and backbone"))
        ref=traj.xyz[0,lig_sel]
        lrmsd=np.sqrt(((traj.xyz[:,lig_sel]-ref)**2).sum(-1).mean(-1))*10
        np.savetxt(f"{out}/ligand_rmsd.csv", np.c_[traj.time/1000, lrmsd], delimiter=",",
                   header="time_ns,ligand_rmsd_A", comments="")
        figs.append(("ligand_rmsd",traj.time/1000,lrmsd,"Ligand RMSD (Å)"))

    # plot
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    n=len(figs)+1
    fig,axes=plt.subplots(1,n,figsize=(4*n,3.2))
    for ax,(nm,x,y,lab) in zip(axes[:len(figs)],figs):
        ax.plot(x,y,lw=0.8); ax.set_xlabel("time (ns)"); ax.set_ylabel(lab)
    axes[-1].plot(resids,rmsf,lw=0.8,color="#c0392b"); axes[-1].set_xlabel("residue"); axes[-1].set_ylabel("RMSF (Å)")
    fig.tight_layout(); fig.savefig(f"{out}/md_analysis.png",dpi=150)
    print(f"analysis -> {out}  (mean bb RMSD {rmsd.mean():.2f} A"
          + (f", mean ligand RMSD {lrmsd.mean():.2f} A" if len(lig_sel)>0 else "") + ")")

if __name__=="__main__":
    main()
'''
open("work/cluster_packages/03_md/scripts/analyze_traj.py", "w").write(analyze_traj_py)

# Write README.md
readme_md = '''# MD simulation package — cluster-ready

Parameterized, solvated, energy-minimized and **equilibrated** protein-ligand
systems for the RORγt lead, plus GPU production + analysis scripts. Built with
**OpenMM 8.5** — Amber **ff14SB** (protein) + **OpenFF-2.1.0** (ligand, via
openmmforcefields) + **TIP3P** water, 0.15 M NaCl, 1.0 nm padding.

## Contents