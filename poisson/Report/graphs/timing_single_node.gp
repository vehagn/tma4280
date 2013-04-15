set term tikz color solid size 5in,3in
unset logscale;
set xrange [0:13];
set xlabel '$P = T_M$'
set ylabel '$S_P$'
set ytics 1 nomirror
set yrange [0.5:12.5]
set y2label '$\eta_P$'
set y2tics 0.1 nomirror;
set y2range [0.5:1]
set output "timing_single_node.tex"
set key right bottom

plot \
"timing_single_node.dat" using 1:2 title '$S_P$' with linespoints linetype 1,\
x title 'Ideal speedup',\
"timing_single_node.dat" using 1:3 title '$\eta_P$' with linespoints axes x1y2;

