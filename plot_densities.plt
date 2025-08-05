set grid
set title "Densities"

set xlabel "mt"
set ylabel "Кин"
set y2label "Пот"

set datafile separator comma

plot "data.csv" using 1:2 with points pt 7 ps 0.7 lc rgb "green"  axes x1y1 title "Kin", "data.csv" using 1:3 with points pt 7 ps 0.7 lc rgb "pink" axes x1y2 title "Pot"

set terminal pngcairo size 1000,1000 enhanced font "Arial,12"
set output "densities.png"
replot
set output

