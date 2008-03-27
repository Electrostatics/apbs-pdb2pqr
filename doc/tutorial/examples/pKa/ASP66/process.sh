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
echo "Calculating the ASH transfer free energy..."
G1=`grep Global 2LZT-ASH66.out | awk '{print $6}'`
G2=`grep Global 2LZT-noASH66.out | awk '{print $6}'`
G3=`grep Global ASH66.out | awk '{print $6}'`
dG_xfer_ASH66=`python -c "print ($G1 - $G2 - $G3)"`
echo "ASH66 transfer free energy FROM TOTAL ENERGY = $dG_xfer_ASH66"

echo "##########################################"
echo "Calculating the ASP transfer free energy..."
G1=`grep Global 2LZT-ASP66.out | awk '{print $6}'`
G2=`grep Global 2LZT-noASP66.out | awk '{print $6}'`
G3=`grep Global ASP66.out | awk '{print $6}'`
dG_xfer_ASP66=`python -c "print ($G1 - $G2 - $G3)"`
echo "ASP66 transfer free energy FROM TOTAL ENERGY = $dG_xfer_ASP66"

echo "##########################################"
echo "Calculating the transfer free energy difference..."
ddG_xfer=`python -c "print ($dG_xfer_ASP66 - $dG_xfer_ASH66)"`
echo "Transfer free energy difference FROM TOTAL ENERGY = $ddG_xfer"

echo "##########################################"
echo "Calculating the pKa shift..."
dpKa=`python -c "print (${ddG_xfer}/${kT}/2.303)"`
echo "pKa shift FROM TOTAL ENERGY = $dpKa"

