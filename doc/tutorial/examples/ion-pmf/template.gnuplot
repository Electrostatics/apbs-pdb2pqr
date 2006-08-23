set size ratio 0.5
set xlabel "Separation (A)"
set ylabel "Total interaction energy (kJ/mol)"
set xrange[-1:11]
set yrange[-30:10]
plot -17.686/x title "Coulomb (dielectric 78.54)" 
set samples 1000
replot 6.6*sqrt(4-(x-0)*(x-0))-31 title 'Sphere 1' with lines
replot 6.6*sqrt(4-(x-SPHERE)*(x-SPHERE))-31 title 'Sphere 2' with lines
replot 6.6*sqrt((1.4*1.4)-(x-WATER)*(x-WATER))-31 title 'Water' with lines
replot 'tot-energy.dat' notitle with linespoints 
