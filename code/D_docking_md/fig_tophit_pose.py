# Auto-extracted generating script
# Produces: fig_tophit_pose.png
# Conda env: dock-md   (run with this environment activated)
# Inputs (expected alongside / in data/): RORC_apo.pdb, CHEMBL3314024.pdb
# Source artifact version: 804e1184-5c52-437f-bb77-39dfab761b56
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import matplotlib.pyplot as plt, matplotlib as mpl, numpy as np
mpl.rcParams.update({"font.size":8,"axes.linewidth":0.6,"font.family":"DejaVu Sans"})

def read_atoms(path, rectype=("ATOM",)):
    ca=[]; allp=[]
    for ln in open(path):
        if ln.startswith(rectype):
            xyz=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
            allp.append(xyz)
            if ln[12:16].strip()=="CA": ca.append(xyz)
    return np.array(ca), np.array(allp)

ca,_=read_atoms("RORC_apo.pdb")
lig=[]
for ln in open("CHEMBL3314024.pdb"):
    if ln.startswith(("ATOM","HETATM")):
        lig.append((float(ln[30:38]),float(ln[38:46]),float(ln[46:54])))
lig=np.array(lig)
# pocket residues: CA within 10A of ligand centroid
ctr=lig.mean(0)
near=ca[np.linalg.norm(ca-ctr,axis=1)<12]

fig=plt.figure(figsize=(7,6))
ax=fig.add_subplot(111,projection='3d')
# backbone
for i in range(len(ca)-1):
    ax.plot(ca[i:i+2,0],ca[i:i+2,1],ca[i:i+2,2],color='0.8',lw=0.6)
# pocket-lining CA
ax.scatter(near[:,0],near[:,1],near[:,2],s=18,c='#3498db',alpha=0.5,edgecolors='none',label='pocket-lining residues')
# ligand
ax.scatter(lig[:,0],lig[:,1],lig[:,2],s=45,c='#e74c3c',edgecolors='k',lw=0.4,zorder=10,label='top hit (CHEMBL3314024)')
# zoom to pocket
r=14
ax.set_xlim(ctr[0]-r,ctr[0]+r);ax.set_ylim(ctr[1]-r,ctr[1]+r);ax.set_zlim(ctr[2]-r,ctr[2]+r)
ax.set_axis_off()
ax.set_title("Top docked hit in RORγt LBD pocket\nCHEMBL3314024, Vina affinity −11.33 kcal/mol",fontsize=8.5)
ax.legend(frameon=False,fontsize=7,loc='upper left')
ax.view_init(elev=15,azim=60)
fig.savefig("fig_tophit_pose.png",dpi=190,bbox_inches='tight')
print("saved top-hit pose figure; ligand centroid", ctr.round(1), "pocket residues", len(near))