#!/bin/bash

CD_LINE
source /usr/software/sourceall.sh
micromamba activate /usr/software/my_env

TOTAL_CPUS=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK)) # Count up the jobs
export SLURM_NTASKS=$TOTAL_CPUS
export CONFIG="docking_config.yaml"
echo "MPI tasks: ${SLURM_NTASKS}"
echo "Config: ${CONFIG}"
python3 /usr/software/dockflow/docking_module/dock6_pipeline/dock6_pipeline.py --config "${CONFIG}" --overwrite
