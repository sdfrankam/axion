set grid
set title "Full density"

set xlabel "mt"
set ylabel "density"

set datafile separator comma

plot "full_energy.csv" using 1:2 with points pt 7 ps 0.7 lc rgb "green"  axes x1y1 title "Full_energy_density"

set terminal pngcairo size 1000,1000 enhanced font "Arial,12"
set output "full_energy.png"
replot
set output

