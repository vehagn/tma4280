set term tikz color solid size 5in,3in
set logscale xy
unset label
unset xrange
set yrange [1.0E-10:3.0E-05]
set format y "%G"
set xlabel '$n$'
set ylabel '$\tau/n^2\log(n)$'
set output "time1.tex"
set key spacing 2

plot \
"time1.dat" using 1:2 title '$N=1$, $M=1$, $T_M=6$, $P=6$' with linespoints, \
"time1.dat" using 1:3 title '$N=1$, $M=2$, $T_M=6$, $P=12$' with linespoints, \
"time1.dat" using 1:4 title '$N=3$, $M=2$, $T_M=6$, $P=36$' with linespoints
