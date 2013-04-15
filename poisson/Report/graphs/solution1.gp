set term tikz color solid size 6in,4in
set output "solution1.tex"

#set terminal png
#set output "heatmap.png"

set cbrange [-1:1];
set xrange [0:63];
set yrange [0:63];
set zrange [-1:1.01];

unset ztics
unset xtics
unset ytics

set ticslevel 0

set samples 2500
set isosamples 2000

set palette defined ( 0 "green", 1 "blue", 2 "red", 3 "orange" )
set view 55,63

set pm3d at s hidden3d 100
set style line 100 linecolor rgb "black"
unset hidden3d
unset surf
unset colorbox
set cbrange [-1:1];

set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )

splot "solution1.dat" matrix with pm3d title ''

