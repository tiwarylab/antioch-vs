#!/bin/bash
#SBATCH --job-name=CTMD
#SBATCH --output=debugMetaD.out # Output File
#SBATCH --error=debugMetaD.err # Error File
#SBATCH --partition=gpu # Need 1 GPU
#SBATCH --gpus=a100:1 # Need 1 GPU
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=8000 # MB
#SBATCH --account=ACCOUNT

