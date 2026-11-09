#!/bin/bash
export AMBERHOME=/programs/x86_64-linux/system/sbgrid_bin
export CUDA_VISIBLE_DEVICES="0"

# SYSTEM
export prmtop=System_Prep/Combine/system.prmtop
export name=System_Prep/Combine/system
export hmass=system_mol1.hmass.prmtop

N_METAD=3

for metadno in `seq 1 $N_METAD`
do
	echo "Running iteration $metadno of $N_METAD"
## 100ns NPT run ##
if test -f "07_Prod.rst"
then
	echo "Production already exists. It will be REMOVED!"
	rm -f 07_Prod.rst 07_Prod.out
fi
pmemd.cuda \
-O \
-i 07_Prod.in \
-o 07_Prod.out \
-p ../$prmtop \
-c 06_Calpha.rst \
-r 07_Prod.rst \
-x 07_Prod.nc \
-inf 07_Prod.mdinfo
if test $? -ne 0
then
	echo "Metadynamics seems to have crashed!"
	echo "This may be intended (such as automatically ending with ligand leaves the protein)"
fi
mv HILLS HILLS_"$metadno"
mv COLVAR COLVAR_"$metadno"
mv 07_Prod.rst 07_Prod_iter"$metadno".rst
mv 07_Prod.out 07_Prod_iter"$metadno".out
mv 07_Prod.nc 07_Prod_iter"$metadno".nc
done
