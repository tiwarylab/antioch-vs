# Dock6 Test Files

This folder contains a small Dock6 flexible-docking test set.

- [`1xgj_protein_nH_ncH.pdb`](1xgj_protein_nH_ncH.pdb) is the input receptor PDB.
- [`1xgj_ligand_nH_ncH.mol2`](1xgj_ligand_nH_ncH.mol2) is the reference ligand used to specify the docking site.
- [`1xgj_min_scored.mol2`](1xgj_min_scored.mol2) is the minimized ligand that should be used for the actual Dock6 docking run.
- [`example_chimera_config.yaml`](example_chimera_config.yaml) is the receptor-preparation config.
- [`example_dock6_config.yaml`](example_dock6_config.yaml) is the Dock6 flexible-docking config.

The DOCK6 pipeline can use either a native installation or a Singularity `.sif`
image. For container execution, set `environment.sif_image` and use
`environment.dock6_dir` for the DOCK6 installation directory inside the image.
The standard internal commands (`sphgen`, `sphere_selector`, `showbox`, `grid`,
and `dock6`) are used automatically. Override `environment.container_commands`
only if your image exposes different command names or paths.
On an HPC cluster, you may first need to run:

```bash
module load singularity
```

Add any extra host directories required by the container under
`environment.singularity_binds`. Input and output directories are bound
automatically.

After docking finishes, load the input receptor PDB, the reference ligand, and
the docking output in your molecular viewer. Comparing those structures helps
you see how well Dock6 flexible docking reproduced the expected binding pose.

For larger screens, tune Dock6 sampling based on your goal. The main docking
module README includes three suggested presets: fast HTVS, balanced HTVS, and
thorough redocking/debugging. In general, higher `max_orientations`, more
preclustered conformers, and final simplex minimization improve sampling at the
cost of runtime.
