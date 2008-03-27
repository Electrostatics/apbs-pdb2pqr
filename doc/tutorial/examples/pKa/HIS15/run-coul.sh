#!/bin/bash 

for i in *.pqr; do
	echo "Running ${i}..."
	coulomb $i | tee ${i%.pqr}-coul.out
done
