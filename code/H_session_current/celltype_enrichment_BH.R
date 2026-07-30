# Auto-extracted generating script
# Produces: celltype_enrichment_BH.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 29ba858c-3b9b-4d55-9ce6-06205c8bf6a6
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

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


import base64
import subprocess

c = host.compute.create("ssh:n003")
base = "/net/dali/home/mscbio/sba50/Soahum/Project/psoriasis/psoriasis-1-bulk"
rscript = r'''
suppressMessages(library(Seurat))
rf <- "RF"
so <- readRDS(file.path(rf,"reference_scissor_full.rds"))
md <- so@meta.data; cl <- md$scissor; ct <- md$celltype
cts <- sort(unique(ct))
res <- lapply(cts, function(k){
  inpos <- cl=="Scissor+"; inbg <- cl=="Background"
  a<-sum(inpos&ct==k); b<-sum(inpos&ct!=k); cc<-sum(inbg&ct==k); d<-sum(inbg&ct!=k)
  ft<-fisher.test(matrix(c(a,b,cc,d),2))
  data.frame(celltype=k, OR=round(ft$estimate,2), p=ft$p.value)
})
res<-do.call(rbind,res)
res$fdr_BH <- p.adjust(res$p, method="BH")
res<-res[order(-res$OR),]
res$p<-signif(res$p,2); res$fdr_BH<-signif(res$fdr_BH,2)
print(res, row.names=FALSE)
write.csv(res, file.path(rf,"celltype_enrichment_BH.csv"), row.names=FALSE)
'''.replace("RF", base+"/results_full")
b64=base64.b64encode(rscript.encode()).decode()
out=c.call_command(
    f"cd {base} && printf %s {b64} | base64 -d > /tmp/enr_bh.R && timeout 55 "
    f"/net/dali/home/mscbio/sba50/miniconda3/envs/scissor-r/bin/Rscript /tmp/enr_bh.R",
    intent="BH-FDR for full-census celltype enrichment")
print(out.get("stdout","")[:1500]); print("ERR:", out.get("stderr","")[-300:])

out=c.call_command(f"cd {base} && base64 results_full/celltype_enrichment_BH.csv", intent="fetch celltype_enrichment_BH.csv")
data=base64.b64decode(out.get("stdout","").replace("\n",""))
open("celltype_enrichment_BH.csv","wb").write(data)
print("celltype_enrichment_BH.csv", len(data), "bytes")