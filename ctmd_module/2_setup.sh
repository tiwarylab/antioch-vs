#!/bin/bash
#

#Set this to "True" to add ions manually
ADD_IONS=$1
parallel -j 12 'echo {} && cd `echo {}|sed s/".mol2"//g` && ../_setup.sh ligand_raw.mol2 `cat Qnet` '$ADD_IONS' && cd ..' ::: `awk -e '{print $1}' hit_data_mol2.spc`
#do
#	echo $f
#	folname=`echo $f|sed s/".mol2"//g`
#	cd $folname
#	../_setup.sh ligand_raw.mol2 `cat Qnet` $ADD_IONS
#	cd ..
#done
