#!/bin/bash

# Put the APBS executeable dir in bindir, i.e:
#bindir="../../../dist/bin/i686-pc-linux/"
bindir=''

if [ "$bindir" = '' ]; then
    echo "Please edit runme.sh to add your APBS executeable directory to the path."
    exit
fi


coord='-2.000 -1.00 0.000 1.000 2.000'

for i in $coord; do
  echo "Moving ion to $i..."
# TO TRANSLATE IN THE X DIRECTION:
  echo "ATOM      1  I   ION     1      -3.000   0.000   0.000  1.00  0.00" \
     > mol1.pqr
  echo "ATOM      1  I   ION     1      $i    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  ${bindir}/apbs apbs.in > OUTPUT_${i} 2>&1
done
