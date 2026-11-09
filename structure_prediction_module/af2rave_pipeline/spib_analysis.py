#!/usr/bin/env python3
# coding: utf-8
#
# spib_analysis.py — env-driven adaptation of prediction/SPIBAnalysis.py.
# Same analysis logic; the hardcoded path roots are replaced by environment variables so run_spib_stage.py
# can drive it, and the AF2 static-structure overlay is guarded so the MD-side analysis runs before those
# static-structure COLVARs exist.
#
# Needs an environment with af2rave, mdtraj, scikit-learn, matplotlib. Driven by run_spib_stage.py; can
# also be run directly from the MD-output dir.
#
# Env vars:
#   SPIB_BASE            folder with <system>/COLVAR subdirs + spib_results/   (default: spib_stage)
#   SPIB_OUT             where protein_outputs_* go                            (default: $SPIB_BASE)
#   SPIB_NUM_FINAL       number of states selected for docking                 (default: 5)
#   AF2_STRUCTURE_ROOT   AF2 static-structure dir (subdirs with *.pdb + *_COLVAR); empty = skip
#                        the AF2-overlay section (FES still produced from MD).

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.covariance import MinCovDet
import mdtraj
import shutil
import tqdm

# ### Read all pickle files and pick optimal dt

BASE_FOLDER = os.environ.get("SPIB_BASE", "spib_stage")
SPIB_FOLDER = os.path.join(BASE_FOLDER, "spib_results")
OUT_ROOT    = os.environ.get("SPIB_OUT", BASE_FOLDER)
GRAPH_DIR   = os.path.join(BASE_FOLDER, "graph_outputs")

all_spib_pickles=sorted([f for f in os.listdir(SPIB_FOLDER) if f.endswith("pkl")])
state_counts=dict()
for pklf in all_spib_pickles:
    model=pickle.load(open(SPIB_FOLDER+"/"+pklf,"rb"))
    model_dt=int(pklf.replace("spib_dt_","").split(".")[0])
    state_counts[model_dt]=model.n_converged_states
time_lags=list(state_counts.keys())
num_states=list(state_counts.values())
order=np.argsort(time_lags)
time_lags=np.array(time_lags)[order]
num_states=np.array(num_states)[order]


optimal_time_lag=time_lags[0]
lowest_state_count=num_states[0]
for i in range(1,len(time_lags)):
    if num_states[i]<lowest_state_count:
        if (lowest_state_count-num_states[i])>(num_states[0]//10):
            optimal_time_lag=time_lags[i]
            break
        lowest_state_count=num_states[i]
print("Optimal time lag:",optimal_time_lag)

if optimal_time_lag==time_lags[0]:
    print("[WARN] Significant state mixing may not have been observed. Consider re-running the Molecular Dynamics part longer")
    optimal_time_lag=time_lags[np.argmin(num_states)].item()
#CV_dt=int(1000)
CV_dt=int(optimal_time_lag)


fig,ax=plt.subplots(figsize=(12,8))
grph=ax.plot(time_lags,num_states,marker='X',markersize=12, color="black",zorder=1)
grph[0].set_markeredgecolor("black")
grph[0].set_markerfacecolor("yellow")
plt.xscale("log")
plt.xticks(time_lags,time_lags,rotation=45,fontsize=16)
plt.yticks(fontsize=16)
plt.xticks(fontfamily='Noto Serif')
plt.yticks(fontfamily='Noto Serif')
ax.set_xlabel("Time lag in frames (logscale)",fontsize=21)
ax.set_ylabel("Number of final states", fontsize=21)
ax.axvline(x=CV_dt,linestyle="--",linewidth=3,alpha=0.7,color="#991111",zorder=0)

os.makedirs(GRAPH_DIR, exist_ok=True)

plt.tight_layout()
plt.savefig(GRAPH_DIR+"/time_lag_analysis.svg")
plt.savefig(GRAPH_DIR+"/time_lag_analysis.png", dpi=300)


# ### Project all COLVARs with this dt

from pathlib import Path
from af2rave import Colvar
import af2rave.spib as af2spib


projection_model=None
def discover_colvars(base: Path) -> list[Path]:
    """Find all COLVAR files in immediate subdirectories."""
    colvars = sorted(base.glob("*/COLVAR"))
    return [c for c in colvars if c.is_file()]
def project_and_save(model_path: str, colvar_path: Path,
                     output_path: Path | None = None, return_directly: bool = False) -> Path | np.ndarray:
    """Load a SPIB model, project a COLVAR, and save the result.

    The output file contains: time, CV1, CV2, ..., CVk, state
    where 'state' is the predicted state number (0-indexed).
    NaN is written for points outside the convex hull of latent space.
    """
    global projection_model
    if projection_model is None:
        projection_model=af2spib.SPIBResult.from_file(model_path)
    result = projection_model

    dt = result.dt
    n_states = int(result.n_converged_states)

    cv = Colvar.from_file(str(colvar_path))

    projected = result.project_colvar(cv)
    n_cvs = projected.shape[0]
    n_frames = projected.shape[1]

    latent_X = projected.data  # (n_cvs, n_frames)
    state_labels = result.project_state_label(latent_X)  # (n_frames,)

    cv_headers = [f"CV{i + 1}" for i in range(n_cvs)]
    all_headers = ["time"] + cv_headers + ["state"]

    if output_path is None:
        out_name = f"spib_cv_dt{dt}.dat"
        output_path = colvar_path.parent / out_name

    if return_directly:
        ret_data=[]
        for frame_idx in range(n_frames):
            cv_vals = projected.data[:, frame_idx]
            state = float(state_labels[frame_idx])
            ret_data.append(np.array(list(cv_vals)+[int(state) if not np.isnan(state) else np.nan]))
        return np.stack(ret_data,axis=0)
    with open(output_path, "w") as f:
        f.write("#! FIELDS " + " ".join(all_headers) + "\n")
        time_arr = projected.time.flatten() if projected.time is not None else np.arange(n_frames)
        for frame_idx in range(n_frames):
            t = float(time_arr[frame_idx])
            cv_vals = projected.data[:, frame_idx]
            state = float(state_labels[frame_idx])
            parts = [f"{t:.6f}"] + [f"{float(v):.6f}" for v in cv_vals]
            if np.isnan(state):
                parts.append("nan")
            else:
                parts.append(f"{int(state)}")
            f.write(" ".join(parts) + "\n")

    print(f"  Saved: {output_path.name}  ({n_cvs} CVs + state x {n_frames} frames)")
    return output_path


MD_FOLDERS=sorted([f for f in os.listdir(BASE_FOLDER) if os.path.isdir(BASE_FOLDER+"/"+f) and ("COLVAR" in list(os.listdir(BASE_FOLDER+"/"+f)))])
assert ("spib_dt_"+str(CV_dt)+".pkl" in all_spib_pickles), "Pickle file for the SPIB model with dt="+str(CV_dt)+" not found"
model_path=SPIB_FOLDER+"/"+"spib_dt_"+str(CV_dt)+".pkl"
for fol in MD_FOLDERS:
    project_and_save(model_path, Path(BASE_FOLDER+"/"+fol+"/COLVAR"))


# ### After dt is picked

trajectories=[]
state_labels=[]
trajectory_lengths=[]
reliable_frames=[]
nan_points=0
total_points=0
for fol in MD_FOLDERS:
    raw_data=np.loadtxt(f"{BASE_FOLDER}/{fol}/spib_cv_dt{CV_dt}.dat")
    reliable_data=(~np.isnan(raw_data[:,-1]))
    trajectories.append(raw_data[reliable_data][:,(1,2)])
    state_labels.append(np.round(raw_data[reliable_data,-1]).astype(int))
    nan_points+=len(raw_data)-np.sum(reliable_data)
    total_points+=len(raw_data)
    # full_data is built from the NaN-dropped rows, so trajectory_lengths must count those same
    # rows, and reliable_frames records their ORIGINAL frame indices. A selected index (into
    # full_data) then maps back to the exact trajectory frame in extract_frame() below. Using the
    # full length here instead offsets every pick after a NaN frame by the cumulative NaN count.
    trajectory_lengths.append(int(np.sum(reliable_data)))
    reliable_frames.append(np.where(reliable_data)[0])

if (nan_points/total_points)>0.01: # More than 1% nan
    print("[WARN]: Many points fall outside canonical state boundaries of SPIB. You might want to run the Molecular Dynamics longer, or run SPIB with higher IB beta")
full_data=np.concatenate(trajectories,axis=0)
full_state_labels=np.concatenate(state_labels,axis=0)


fig,ax=plt.subplots(1,2,figsize=(12,5))

for ti,t in enumerate(trajectories):
    ax[0].scatter(t[:,0],t[:,1],s=3,color="C"+str(ti))
ax[0].set_xticks([],[])
ax[0].set_yticks([],[])
ax[0].set_xlabel("$\\sigma_0$",fontsize=21)
ax[0].set_ylabel("$\\sigma_1$",fontsize=21)
ax[0].set_title("Trajectory Mixing",fontsize=24)

ax[1].scatter(full_data[:,0],full_data[:,1],s=3,c=full_state_labels)
ax[1].set_xticks([],[])
ax[1].set_yticks([],[])
ax[1].set_xlabel("$\\sigma_0$",fontsize=21)
ax[1].set_ylabel("$\\sigma_1$",fontsize=21)
ax[1].set_title("Converged labels",fontsize=24)

os.makedirs(GRAPH_DIR, exist_ok=True)

plt.tight_layout()
plt.savefig(GRAPH_DIR+"/converged_states.svg")
plt.savefig(GRAPH_DIR+"/converged_states.png", dpi=300)


# ### (Pseudo)-Free Energy Surface

kT = 0.593 # kcal/mol
GRID_SIZE=64
FES_GRID=np.zeros((GRID_SIZE,GRID_SIZE),dtype=float)
x_min=np.min(full_data[:,0])
x_max=np.max(full_data[:,0])
y_min=np.min(full_data[:,1])
y_max=np.max(full_data[:,1])

x_bins=np.linspace(x_min,x_max,GRID_SIZE)
y_bins=np.linspace(y_min,y_max,GRID_SIZE)
BATCH_SIZE=4096

x_record=[]
y_record=[]

start_idx=0
while start_idx<len(full_data):
    subdata=full_data[start_idx:start_idx+BATCH_SIZE]
    start_idx+=BATCH_SIZE
    if len(subdata)<=1: break
    x_idxs=np.digitize(subdata[:,0],bins=x_bins)-1
    y_idxs=np.digitize(subdata[:,1],bins=y_bins)-1
    np.add.at(FES_GRID,(x_idxs,y_idxs),1)
    x_record.append(x_idxs)
    y_record.append(y_idxs)

x_record=np.concatenate(x_record)
y_record=np.concatenate(y_record)

FES_GRID=(FES_GRID[:,1:]+FES_GRID[:,:-1])/2
FES_GRID=(FES_GRID[1:,:]+FES_GRID[:-1,:])/2
x_record_retro=np.round(x_record.astype(float)*((GRID_SIZE-1)/GRID_SIZE)).astype(int)
y_record_retro=np.round(y_record.astype(float)*((GRID_SIZE-1)/GRID_SIZE)).astype(int)

F = -kT * np.log(FES_GRID)
FES_FOCUS=4*kT
FES_RESOLVE=5
F -= F.min()


# ### Identifying stable minimas

NUM_FINAL_STRUCTURES=int(os.environ.get("SPIB_NUM_FINAL", 5))
EXCESS=0
from sklearn.mixture import GaussianMixture
GMM_SEED=int(os.environ.get("SPIB_GMM_SEED","9"))   # pin so the selected docking structures are reproducible
gmix=GaussianMixture(n_components=NUM_FINAL_STRUCTURES+EXCESS, random_state=GMM_SEED)
gaussian_class=gmix.fit_predict(full_data)

plt.figure(figsize=(8,8))
plt.scatter(full_data[:,0],full_data[:,1],c=gaussian_class,s=3)
plt.xlabel("$\\sigma_0$",fontsize=24)
plt.ylabel("$\\sigma_1$",fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title(f"{NUM_FINAL_STRUCTURES+EXCESS} states sampled for docking",fontsize=28)
plt.savefig(GRAPH_DIR+"/docking_selection.svg")
plt.savefig(GRAPH_DIR+"/docking_selection.png", dpi=300)


point_FEs=F[x_record_retro,y_record_retro]
gaussian_FEs=dict()
for ug in np.unique(gaussian_class):
    ener=point_FEs[gaussian_class==ug]
    wts=np.exp(-ener/kT)
    ener_val=np.average(ener,weights=wts)
    gaussian_FEs[ug]=ener_val
gids=np.array(list(gaussian_FEs.keys()))
geners=np.array(list(gaussian_FEs.values()))
gids=gids[np.argsort(geners)]
geners=geners[np.argsort(geners)]

selected_minimas=[]
selected_minimas_idx=[]
for ug in gids:
    reqfilt=((gaussian_class==ug))
    cidx=np.argmin(point_FEs[reqfilt])
    selected_minimas_idx.append(np.where(reqfilt)[0][cidx])
    selected_minimas.append(full_data[selected_minimas_idx[-1]])

selected_minimas=np.array(selected_minimas)
selected_minimas_idx=np.array(selected_minimas_idx)


# ### Project and pick  (AF2 static-structure overlay — optional)

AF2_STRUCTURE_ROOT=os.environ.get("AF2_STRUCTURE_ROOT","")
HAVE_AF2 = bool(AF2_STRUCTURE_ROOT) and os.path.isdir(AF2_STRUCTURE_ROOT)
structure_SPIBs=[]
all_pdb_files=[]
if not HAVE_AF2:
    print("[info] AF2_STRUCTURE_ROOT unset/missing — skipping AF2 static-structure overlay "
          "(set it to a dir of AF2 PDBs + *_COLVAR to enable).")
else:
    for fol in os.listdir(AF2_STRUCTURE_ROOT):
        if not os.path.isdir(AF2_STRUCTURE_ROOT+"/"+fol): continue
        else: print(fol)
        ALL_AF2_STRUCTURES=AF2_STRUCTURE_ROOT+"/"+fol
        pdb_files=sorted([f for f in os.listdir(ALL_AF2_STRUCTURES) if f.endswith("pdb")])
        all_pdb_files+=[ALL_AF2_STRUCTURES+"/"+p for p in pdb_files]
        colvar_files=[f.replace(".pdb","_COLVAR") for f in pdb_files]
        print("\t",len(pdb_files),"AF2 structures found")

        if os.path.exists(ALL_AF2_STRUCTURES+"/projection.npy"):
            structure_SPIBs=np.load(ALL_AF2_STRUCTURES+"/projection.npy")
            print("\t Loading data from",ALL_AF2_STRUCTURES+"/projection.npy")
        else:
            for cvf in tqdm.tqdm(colvar_files,ncols=100):
                proj_data=project_and_save(model_path, Path(ALL_AF2_STRUCTURES+"/"+cvf),Path("./temp_proj.dat"),return_directly=True)
                structure_SPIBs.append(proj_data.squeeze())
            np.save(ALL_AF2_STRUCTURES+"/projection.npy",structure_SPIBs)
    if not len(structure_SPIBs): raise FileNotFoundError(f"No PDB files found in {AF2_STRUCTURE_ROOT}")
    structure_SPIBs=np.stack(structure_SPIBs,axis=0)
    print(structure_SPIBs.shape[0],"total structures collected")


# --- FES plot ---
fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    F.T,
    origin='lower',
    aspect='auto',
    cmap='seismic',
    extent=[x_min, x_max, y_min, y_max],
    vmin=0,
    vmax=FES_FOCUS+2*kT,
    alpha=0.95
)

cs = ax.contour(
    F.T,
    levels=np.arange(F.min(), FES_FOCUS + (FES_FOCUS/FES_RESOLVE), FES_FOCUS/FES_RESOLVE),
    colors='black',
    linewidths=1.25,
    alpha=0.8,
    extent=[x_min, x_max, y_min, y_max],
)

cbar = fig.colorbar(im, ax=ax, pad=0.02)
cbar.set_label('Free Energy (kcal/mol)', fontsize=12)

ax.set_xlabel('$\\sigma_0$', fontsize=16)
ax.set_ylabel('$\\sigma_1$', fontsize=16)
ax.set_xticks([],[])
ax.set_yticks([],[])
ax.set_title('Free Energy Surface', fontsize=18)

ax.scatter(selected_minimas[:NUM_FINAL_STRUCTURES,0],selected_minimas[:NUM_FINAL_STRUCTURES,1],marker='o',color="#DDAA00",label="Selected States", edgecolor="red", linewidth=1)
ax.scatter(selected_minimas[NUM_FINAL_STRUCTURES:,0],selected_minimas[NUM_FINAL_STRUCTURES:,1],marker='o',color="#FF11CC", label="Rejected States")
ax.legend(fontsize=18)

if HAVE_AF2 and len(structure_SPIBs):
    ax.scatter(structure_SPIBs[:,0],structure_SPIBs[:,1],color="cyan",s=10,alpha=0.6)

plt.tight_layout()
plt.savefig(GRAPH_DIR+"/FES.svg")
plt.savefig(GRAPH_DIR+"/FES.png", dpi=300)


# --- copy the nearest AF2 static structure to each selected minimum (only if AF2 overlay ran) ---
if HAVE_AF2 and len(structure_SPIBs):
    direct_out=os.path.join(OUT_ROOT,"protein_outputs_direct")
    os.makedirs(direct_out, exist_ok=True)
    for sm in selected_minimas:
        af2_idx=np.argmin(np.sum((structure_SPIBs[:,:2]-sm[np.newaxis,:])**2,axis=-1))
        shutil.copy(all_pdb_files[af2_idx],direct_out+"/")


# --- extract the representative MD frame for each selected minimum ---
def extract_frame(direct_idx):
    trajidx=-1
    for i in range(len(trajectory_lengths)):
        if direct_idx-trajectory_lengths[i]<0:
            trajidx=i
            break
        direct_idx-=trajectory_lengths[i]

    if trajidx<0: raise ValueError(f"Index {direct_idx} can't be found!")
    else: print("Trajectory number:",trajidx+1)
    # direct_idx is now the position among this trajectory's NaN-dropped rows; map it back to the
    # original frame so the saved structure is exactly the SPIB-selected frame (not NaN-offset).
    orig_frame=int(reliable_frames[trajidx][direct_idx])
    traj_path=BASE_FOLDER+"/"+MD_FOLDERS[trajidx]+"/"
    traj_data=mdtraj.load_xtc(traj_path+"/06p5_Prod.xtc",traj_path+"/only_protein.pdb")
    return traj_data[orig_frame]

md_out=os.path.join(OUT_ROOT,"protein_outputs_MD")
os.makedirs(md_out, exist_ok=True)
for finidx in selected_minimas_idx:
    extract_frame(finidx).save_pdb(f"{md_out}/protein{finidx}.pdb")

print("Done. Outputs in", GRAPH_DIR, "and", md_out,
      ("(+ "+os.path.join(OUT_ROOT,'protein_outputs_direct')+")" if HAVE_AF2 else "(AF2 overlay skipped)"))
