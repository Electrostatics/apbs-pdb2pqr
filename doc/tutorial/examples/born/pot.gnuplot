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
replot 'potential-1.00-xaxis.dat' title "1.00 A radius" with linespoints pointsize 2
replot 'potential-2.00-xaxis.dat' title "2.00 A radius" with linespoints pointsize 2
replot 'potential-4.00-xaxis.dat' title "4.00 A radius" with linespoints pointsize 2
replot 'potential-6.00-xaxis.dat' title "6.00 A radius" with linespoints pointsize 2
replot 'potential-8.00-xaxis.dat' title "8.00 A radius" with linespoints pointsize 2
replot 'potential-10.00-xaxis.dat' title "10.00 A radius" with linespoints pointsize 2
