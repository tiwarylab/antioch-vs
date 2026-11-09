#!/bin/bash
#SBATCH --job-name=ANTIOCH_MD
#SBATCH --output=debugMD.out
#SBATCH --error=debugMD.err
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1
#SBATCH --time=27:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32000


source /Path/To/singularity_images/ambertools24_plmd/sourceme.sh
folname=`basename $(pwd)`
echo $folname
singularity exec --nv --bind ../:/scratch $AMBERTOOLS24_IMAGE_PATH /scratch/$folname/amber_MD_template/singularity_passer.sh $folname 0

# Setup reference for PLUMED
cpptraj -i amber_MD_template/extract.in
vmd -dispdev text -e amber_MD_template/weightref.tcl
cp amber_MD_template/plumed_template.dat plumed_bpmetad.dat # Change this line to change how plumed is set up

# Actual BPMetad
singularity exec --nv --bind ../:/scratch $AMBERTOOLS24_IMAGE_PATH /scratch/$folname/amber_MD_template/singularity_passer.sh $folname 1
