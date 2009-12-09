#!/bin/bash

# Generate the input file
cat > apbs-pmf.in << EOF

#############################################################################
### BORN ION SOLVATION ENERGY
### $Id: apbs.in 998 2006-11-28 21:24:40Z sobolevnrm $
###
### Please see APBS documentation (http://apbs.sourceforge.net/doc/) for 
### input file sytax.
#############################################################################

# READ IN MOLECULES
read 
    mol pqr ion.pqr       # Read molecule 1
end

# CALCULATE POTENTIAL FOR SOLVATED STATE
elec name solvated
   mg-manual
    dime 65 65 65
    grid 0.33 0.33 0.33
    gcent mol 1
    mol 1 
    lpbe
    bcfl mdh
    pdie 1.0
    sdie 78.54
    chgm spl2
    srfm mol
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    gamma 0.105
    calcenergy total
    calcforce no
end

# CALCULATE POTENTIAL FOR REFERENCE STATE
elec name reference
    mg-manual
    dime 65 65 65
    grid 0.33 0.33 0.33
    gcent mol 1
    mol 1
    lpbe
    bcfl mdh
    pdie 1.0
    sdie 1.0
    chgm spl2
    srfm mol
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    gamma 0.105
    calcenergy total
    calcforce no
end

# COMBINE TO FORM SOLVATION ENERGY
print elecEnergy solvated - reference end
#print apolForce solvated - asdf2 end
#print apolEnergy solvated + asdf2 end

# SO LONG
quit
EOF

apbspath=`which apbs`
if test -z "${apbspath}"; then
	echo ""
	echo "Please add the apbs binary to your path..."
	echo ""
	exit
fi

radii='1.00 2.00 3.00 4.00 5.00 6.00'

for i in $radii; do
  echo "Generating ion with $i radius..."
  echo "ATOM      1  I   ION     1       0.000   0.000  0.000  1.00  $i" \
    > ion.pqr
  apbs apbs-pmf.in > OUTPUT-${i} 2>&1 
done
