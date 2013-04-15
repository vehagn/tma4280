set term tikz color solid size 5in,3in
unset logscale;
set xrange [0:13];
set xlabel '$M$'
set ylabel '$\eta_P$'
set ytics 0.1 nomirror
set output "timing_MPI.tex"
set key right bottom

plot \
"timing_MPI.dat" using 1:6 title '$S_P$' with linespoints;


