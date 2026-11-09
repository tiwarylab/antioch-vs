#!/bin/bash
#SBATCH --job-name=<job_name>
#SBATCH --partition=<partition>
#SBATCH --qos=<qos>
#SBATCH --time=24:00:00
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1
#SBATCH --output=/path/to/logs/decoys_slurm_%j.log
#SBATCH --error=/path/to/logs/decoys_slurm_%j.err
