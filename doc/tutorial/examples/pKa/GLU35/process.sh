#!/bin/bash

# Assumes APBS 0.5.0

epsp=20
kT=2.5

echo ""
echo "##########################################"
echo "# pKa shift from TOTAL ENERGY            #"
echo "##########################################"
echo ""
echo "##########################################"
echo "Calculating the GLH transfer free energy..."
G1=`grep Global 2LZT-GLH35.out | awk '{print $6}'`
G2=`grep Global 2LZT-noGLH35.out | awk '{print $6}'`
G3=`grep Global GLH35.out | awk '{print $6}'`
dG_xfer_GLH35=`python -c "print ($G1 - $G2 - $G3)"`
echo "GLH35 transfer free energy FROM TOTAL ENERGY = $dG_xfer_GLH35"

echo "##########################################"
echo "Calculating the GLU transfer free energy..."
G1=`grep Global 2LZT-GLU35.out | awk '{print $6}'`
G2=`grep Global 2LZT-noGLU35.out | awk '{print $6}'`
G3=`grep Global GLU35.out | awk '{print $6}'`
dG_xfer_GLU35=`python -c "print ($G1 - $G2 - $G3)"`
echo "GLU35 transfer free energy FROM TOTAL ENERGY = $dG_xfer_GLU35"

echo "##########################################"
echo "Calculating the transfer free energy difference..."
ddG_xfer=`python -c "print ($dG_xfer_GLU35 - $dG_xfer_GLH35)"`
echo "Transfer free energy difference FROM TOTAL ENERGY = $ddG_xfer"

echo "##########################################"
echo "Calculating the pKa shift..."
dpKa=`python -c "print (${ddG_xfer}/${kT}/2.303)"`
echo "pKa shift FROM TOTAL ENERGY = $dpKa"

