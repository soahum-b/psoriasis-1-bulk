# Auto-extracted generating script
# Produces: gradient_program_DE_BH.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds
# Source artifact version: 77c061df-01c8-4c00-9518-2d1cf766b3e9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import base64

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")

    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


c = host.compute.create("ssh:n003")
base = "/net/dali/home/mscbio/sba50/Soahum/Project/psoriasis/psoriasis-1-bulk"

rscript = r'''
rf <- "RF"
de <- read.csv(file.path(rf,"gradient_program_DE.csv"))
# protocol correction: BH-FDR over the tested gene family (replaces Seurat Bonferroni p_val_adj)
de$fdr_BH <- p.adjust(de$p_val, method="BH")
write.csv(de, file.path(rf,"gradient_program_DE_BH.csv"), row.names=FALSE)
cat("N_TESTED_FAMILY:", nrow(de), "\n")
cat("SIG_Bonferroni(p_val_adj<0.05):", sum(de$p_val_adj<0.05, na.rm=TRUE), "\n")
cat("SIG_BH(fdr_BH<0.05):", sum(de$fdr_BH<0.05, na.rm=TRUE), "\n")
cat("STAT3_in_family:", "STAT3" %in% de$gene, "\n")
r <- de[de$gene=="STAT3",]; cat("STAT3_ROW:\n"); print(r)
# range check
cat("max p_val among BH-sig:", if(any(de$fdr_BH<0.05)) max(de$p_val[de$fdr_BH<0.05]) else NA, "\n")
'''.replace("RF", base+"/results_full")

b64=base64.b64encode(rscript.encode()).decode()
out = c.call_command(
    f"cd {base} && printf %s {b64} | base64 -d > /tmp/bh2.R && "
    f"timeout 50 /net/dali/home/mscbio/sba50/miniconda3/envs/scissor-r/bin/Rscript /tmp/bh2.R",
    intent="apply BH-FDR to existing gradient DE family and check STAT3")
print(out.get("stdout","")[:1500]); print("ERR:", out.get("stderr","")[-400:])

# fetch result
out=c.call_command(f"cd {base} && base64 results_full/gradient_program_DE_BH.csv | split -b 50000 - /tmp/gde_ && ls /tmp/gde_*",
                   intent="split base64 into chunks")
chunks=out.get("stdout","").split()
b64=""
for ch in chunks:
    o=c.call_command(f"cat {ch}", intent=f"cat {ch}")
    b64+=o.get("stdout","").replace("\n","")
data=base64.b64decode(b64)
open("gradient_program_DE_BH.csv","wb").write(data)
print("full bytes:", len(data), "rows:", data.decode().count(chr(10)))