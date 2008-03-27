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
echo "Calculating the HSP transfer free energy..."
G1=`grep Global 2LZT-HSP15.out | awk '{print $6}'`
G2=`grep Global 2LZT-noHSP15.out | awk '{print $6}'`
G3=`grep Global HSP15.out | awk '{print $6}'`
dG_xfer_HSP15=`python -c "print ($G1 - $G2 - $G3)"`
echo "HSP15 transfer free energy FROM TOTAL ENERGY = $dG_xfer_HSP15"

echo "##########################################"
echo "Calculating the HIS transfer free energy..."
G1=`grep Global 2LZT-HIS15.out | awk '{print $6}'`
G2=`grep Global 2LZT-noHIS15.out | awk '{print $6}'`
G3=`grep Global HIS15.out | awk '{print $6}'`
dG_xfer_HIS15=`python -c "print ($G1 - $G2 - $G3)"`
echo "HIS15 transfer free energy FROM TOTAL ENERGY = $dG_xfer_HIS15"

echo "##########################################"
echo "Calculating the transfer free energy difference..."
ddG_xfer=`python -c "print ($dG_xfer_HIS15 - $dG_xfer_HSP15)"`
echo "Transfer free energy difference FROM TOTAL ENERGY = $ddG_xfer"

echo "##########################################"
echo "Calculating the pKa shift..."
dpKa=`python -c "print (${ddG_xfer}/${kT}/2.303)"`
echo "pKa shift FROM TOTAL ENERGY = $dpKa"

