set xlabel "Radius (A)"
set ylabel "Charge density (e M)"
plot 'qdens.dat' \
  using 1:3 \
  title 'Cumulative charge density'  \
  with linespoints 
replot 5.49 title '76% level'
