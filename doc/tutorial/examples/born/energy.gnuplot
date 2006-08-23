set xlabel "Ion radius (A)"
set ylabel "Solvation energy (kJ/mol)"
plot 'energy.dat' \
    title "Numerical" \
    with linespoints pointsize 2
replot -685.695/x \
    title "Analytical"
    
