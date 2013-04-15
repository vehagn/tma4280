set term tikz color solid size 5in,3.5in
set logscale x
unset label
unset xrange
set xlabel '$n$'
set ylabel '$\eta_P$'
set output "peff.tex"
set key spacing 2
set key out bot center

plot \
"peff.dat" using 1:2 title '$N=3$, $M=12$, $T_M=1$, $P=36$' with linespoints, \
"peff.dat" using 1:3 title '$N=3$, $M=3$, $T_M=4$, $P=36$' with linespoints
