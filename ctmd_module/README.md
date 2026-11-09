# CTMD Module — AMBER MD Setup Automation

Automates the full AMBER preparation pipeline for ligand-receptor complexes: charge calculation, topology generation, solvation, and neutralization.

## Quick Start

### 1. Prepare your input files

Create one working directory per selected receptor-ligand pair. This matters for
ensemble docking: a consensus ligand may choose a different best receptor than
the next ligand, while CTMD uses one receptor PDB per project.

Each project has this structure:

```
your_project/
├── molecules/
│   └── ligand_pose.mol2
└── rec_prepped.pdb          # Prepared receptor PDB file
```

- **`molecules/`** — Place one or more docked ligand poses in MOL2 format here. Names are unrestricted; `prepare_ctmd.sh` creates its charge list from the files it finds.
- **`rec_prepped.pdb`** — Your prepared receptor structure (protonated, missing atoms filled, etc.).

### 2. Activate the AMBER environment

```bash
source /scratch/software/micromamba_root/sourceme.sh
# or fallback:
# export MAMBA_ROOT_PREFIX=/scratch/software/micromamba_root
# eval "$(/scratch/software/micromamba_install/bin/micromamba shell hook -s bash)"
# micromamba activate chem_env
```

### 3. Run the preparation script

From the CTMD project directory (not from `ctmd_module/`):

```bash
bash /path/to/antioch-vs/ctmd_module/prepare_ctmd.sh ./rec_prepped.pdb [add_ions]
```

| Argument | Description | Default |
|----------|-------------|---------|
| `rec_prepped.pdb` | Path to receptor PDB (copied to `./rec_prepped.pdb` if different) | Required |
| `add_ions` | Add Na+/Cl- ions after solvation | `False` |

### 4. Check results

Each ligand gets its own folder containing:

| File | Description |
|------|-------------|
| `system.pdb` | Solvated complex coordinates |
| `system.prmtop` | AMBER topology file |
| `system.inpcrd` | AMBER coordinate file |
| `protlig_vac.*` | Gas-phase files (for MMGBSA) |

Individual logs: `<ligand>_setup.log`
Master log: `ctmd_setup.log`

## What the Script Does

1. **Pre-flight checks** — Verifies all required tools (antechamber, parmchk2, tleap, sqm, parmed) are available
2. **Charge computation** — Sums atomic charges from each mol2 file, generates `hit_data_mol2.spc`
3. **Script preparation** — Copies and patches CTMD scripts for local compatibility (mawk, missing `mul`)
4. **Protein charge** — Computes net charge via tleap + parmed, then patches the receptor-neutralization lines in `setup.leap`
5. **Directory setup** — Creates per-ligand folders with MOL2 copies and Qnet files
6. **AMBER pipeline** — Runs antechamber → parmchk2 → tleap for each ligand (up to 12 parallel, 1 CPU each)

## Requirements

- **AMBER chem_env** — antechamber, parmchk2, tleap, sqm, parmed
- **python3** with parmed module
- **~30 min per ligand** for SQM charge optimization (up to 12 ligand jobs in parallel)

## Troubleshooting

| Error | Fix |
|-------|-----|
| "Required tools not found" | Source chem_env before running |
| "tleap failed to generate prmtop" | Check rec_prepped.pdb for bad residues/missing atoms |
| Individual ligand fails | Check `<ligand>_setup.log` for details |
