set xlabel "X-axis position (A)"
set ylabel "Electrostatic potential (kT/e)"
plot [0:11] [0:350] \
    17.686/2.5/abs(x) \
    title "Coulomb (dielectric 78.54)" \
    with linespoints pointsize 2
replot \
    17.686*78.54/2.5/abs(x) \
    title "Coulomb (dielectric 1.00)" \
    with linespoints pointsize 2
