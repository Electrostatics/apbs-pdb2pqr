#!/bin/bash

if test -z "${APBS}"; then
	echo "Please define the APBS environmental variable to point to your APBS"
	echo "binary.  Under bash, you would type:"
	echo "    export APBS=/path/to/binary"
	exit
fi

coord='-3.000 -2.750 -2.500 -2.250 -2.000 -1.750 -1.50 -1.250 -1.00 -0.750 
-0.500 -0.250 0.000 0.250 0.500 0.750 1.000 1.250 1.500 1.750 2.000 
2.250 2.500 2.750 3.000'

for i in $coord; do
	echo "Moving ion to $i..."
	# TO TRANSLATE IN THE X DIRECTION:
	echo "ATOM      1  I   ION     1      -3.000   0.000   0.000  1.00  2.00" > mol1.pqr
	echo "ATOM      1  I   ION     1      $i    0.000   0.000  1.00 2.00" > mol2.pqr
	cat mol1.pqr mol2.pqr > complex.pqr
	${APBS} apbs.in > OUTPUT_${i} 2>&1
done
