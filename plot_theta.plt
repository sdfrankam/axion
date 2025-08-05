set grid
set title "График theta(mt)"

set xlabel "mt"
set ylabel "theta"

set datafile separator comma

plot "teta_q.csv" using 2:1  with points pt 7 ps 0.7 lc rgb "blue" 
#lc rgb "blue", -0.125*x+3.0 with lines lw 2 lc rgb "green"  

set terminal pngcairo size 800,600 enhanced font "Arial,12"
set output "theta_q.png"
replot
set output
