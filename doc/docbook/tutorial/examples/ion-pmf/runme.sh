#!/bin/bash

rm -f solv-energy.dat
rm -f coul-energy.dat
rm -f tot-energy.dat
rm -f apbs-*.out
rm -f apbs-*.err
rm -f plot-*.gnuplot
rm -f *.dat

for dist in 9.00 8.00 7.00 6.67 6.33 6.00 5.00 4.75 4.50 4.25 4.00 3.75 3.50 3.00 2.50; do

    infile=apbs.in
    errfile=apbs-${dist}.err
    outfile=apbs-${dist}.out
    potfile=potential-${dist}-xaxis.dat

    # Set up structures
    echo "*********************************"
    echo "Separation = ${dist} A"
    echo "$dist" | awk '{printf("ATOM 1 I ION 1 0.000 0.000 0.000 -1.000 2.000\n")}' > mol1.pqr
    echo "$dist" | awk -v d=$dist '{printf("ATOM 1 I ION 1 %4.3f 0.000 0.000 1.000 2.000\n", d)}' > mol2.pqr
    cat mol1.pqr mol2.pqr > complex.pqr

    # Run APBS
    echo "Running APBS..."
    apbs ${infile} > ${outfile} 2> ${errfile}

    # Project potentials
    echo "Projecting potentials..."
    python potential.py pot-solv.dx > pot-solv-xaxis-${dist}.dat
    python potential.py pot-vac.dx > pot-solv-vac-${dist}.dat
    paste pot-solv-xaxis-${dist}.dat pot-solv-vac-${dist}.dat | awk '{print $1, ($4-$2)}' > rxnfield-xaxis-${dist}.dat

    # Get energy
    dGsolv=`grep "Global net energy" ${outfile} | awk '{print $5}'`
    echo "Solvation energy = ${dGsolv} kJ/mol"
    echo ${dist} ${dGsolv} >> solv-energy.dat
    dGcoul=`echo "17.686 ${dist}" | awk '{print -$1*78.54/$2}'`
    echo "Coulomb energy = ${dGcoul} kJ/mol"
    echo ${dist} ${dGcoul} >> coul-energy.dat
    dGtot=`echo ${dGsolv} ${dGcoul} | awk '{print $1+$2}'`
    echo "Total energy = ${dGtot} kJ/mol"
    echo ${dist} ${dGtot} >> tot-energy.dat

    killall -TERM gnuplot_x11 2>&1 > /dev/null
    wdist=`echo ${dist} | awk '{print $1/2}'`
    cat template.gnuplot | \
      sed -e "s/SPHERE/${dist}/g" | \
      sed -e "s/WATER/${wdist}/g" > plot-${dist}.gnuplot
    gnuplot -persist plot-${dist}.gnuplot 2>&1 > /dev/null & 

done


