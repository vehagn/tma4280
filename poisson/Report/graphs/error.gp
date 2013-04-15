set term tikz color solid size 5in,3in
set logscale xy
set xlabel '$n$'
set ylabel 'Error'
set xrange [10:1800]
set format y "%G"
set output "error.tex"
set key spacing 2

plot \
"error.dat" using 1:2 title '$N=1$, $M=1$, $T_M=6$, $P=6$' with linespoints, \
"error.dat" using 1:3 title '$N=3$, $M=2$, $T_M=6$, $P=36$' with linespoints, \
x**-2 title '$n^{-2}$'

