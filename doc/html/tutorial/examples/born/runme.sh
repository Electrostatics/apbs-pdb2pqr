#!/bin/bash 

rm -f apbs-*.in
rm -f apbs-*.out
rm -f apbs-*.err
rm -f ion-*.pqr
rm -f energy.dat

cat coulomb.gnuplot > pot.gnuplot
gnuplot -persist pot.gnuplot 2>&1 > /dev/null &

# Run a series of calculations
for radius in 1.00 2.00 4.00 6.00 8.00 10.00; do

    infile=apbs-${radius}.in
    errfile=apbs-${radius}.err
    outfile=apbs-${radius}.out
    pqrfile=ion-${radius}.pqr
    engfile=energy.dat
    dxfile=potential-${radius}.dx
    potfile=potential-${radius}-xaxis.dat

    echo "*********************************"
    echo "Born radius = ${radius} A" 

    # Set up structures
    cat example-input.in | sed -e "s/RADIUS/${radius}/g" > ${infile}
    cat template.pqr | sed -e "s/RADIUS/${radius}/g" > ${pqrfile}

    # Run APBS
    echo "Running APBS..."
    apbs ${infile} > ${outfile} 2> ${errfile}

    # Get energy
    dG=`grep Global ${outfile} | awk '{printf("%.3f\n", $(NF-1))}'`
    echo ${radius} ${dG} >> ${engfile}
    echo "Solvation energy = ${dG} kJ/mol"

    # Plot it
    echo "Writing x-axis potential values to ${potfile}..."
    python2 potential.py ${dxfile} > ${potfile}
    echo "replot '${potfile}' title \"${radius} A radius\" with linespoints pointsize 2" >> pot.gnuplot
    killall -TERM gnuplot_x11 2>&1 > /dev/null
    gnuplot -persist pot.gnuplot 2>&1 > /dev/null &
done
echo "*********************************"
echo "Plotting solvation energies"
gnuplot -persist energy.gnuplot 2>&1 > /dev/null &
