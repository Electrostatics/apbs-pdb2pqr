#!/bin/bash

apbspath=`which apbs`
if test -z "${apbspath}"; then
	echo ""
	echo "Please add the apbs binary to your path..."
	echo ""
	exit
fi


alkanes="2-methylbutane.pdb butane.pdb cyclohexane.pdb cyclopentane.pdb ethane.pdb hexane.pdb isobutane.pdb methane.pdb neopentane.pdb pentane.pdb propane.pdb"

for pdb in $alkanes; do
	alkane=${pdb%.pdb}
	echo "Calculating for ${alkane}..."
	cp ${pdb} mol.pdb
	apbs apbs-apolar.in > OUTPUT-${alkane} 2>&1 
done
