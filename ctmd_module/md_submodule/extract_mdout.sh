#!/bin/bash

folder=$1
if test -d "$folder"
then
    echo "Folder found: $folder"
else
    echo "Folder $folder does not exist"
    exit 1
fi

DOCKFLOW_ROOT=/usr/software/dockflow
cd $folder
if test -f "06p5_Prod.nc"
then
    cat $DOCKFLOW_ROOT/ctmd_module/md_submodule/strip.in|sed s/PRMTOP_OUT/"system_stripped.prmtop"/|sed s/PRMTOP/"system.prmtop"/|sed s/"NETCDF_TRAJ"/"06p5_Prod.nc"/|sed  s/"XTC_OUT"/"06p5_Prod.xtc"/ > strip_md.in
    cpptraj -i strip_md.in
    cat $DOCKFLOW_ROOT/ctmd_module/md_submodule/framedump.in|sed s/PRMTOP/"system_stripped.prmtop"/|sed s/"XTC_TRAJ"/"06p5_Prod.xtc"/|sed s/"PDB_OUT"/"only_protein.pdb"/ > mdframe.in
    cpptraj -i mdframe.in
else
    echo "In folder $folder: MD trajectory 06p5_Prod.nc was not found"
    exit 1
fi
