set xlabel "Radius (A)"
set ylabel "Charge density (e M)"
plot 'qdens.dat' \
  using 1:2 \
  title 'Charge density' \
  with linespoints
