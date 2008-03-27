#!/bin/bash

for i in *.in; do
	echo "Running ${i}..."
	(apbs $i 2>&1) | tee ${i%.in}.out
done
