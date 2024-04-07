set term postscript enhanced eps "Helvetica" 20
set out "CO_rates.eps"

set xrange [0:0.8]
set xtics 0.1
set mxtics 10
set xlabel "Z/R"

set logscale y
set yrange [0.01:1]
set ytics 10
set mytics 10
set ylabel "Shielding Factor"
set format y "10^{%T}"

set logscale y2
set y2range [1e-14:1e-2]
set y2tics 10
set my2tics 10
set y2label "CO Photodissociation Rate (s^{-1})"
set format y2 "10^{%T}"

set key bottom right

plot 'CO_rates_1AU.dat' u ($1/1):($2)        axis x1y1 w lines lw 3 lt 1 lc 1 t "1 AU",  \
     'CO_rates_1AU.dat' u ($1/1):($3/$2)     axis x1y2 w lines lw 3 lt 3 lc 1 t "",  \
     'CO_rates_1AU.dat' u ($1/1):($3)        axis x1y2 w lines lw 3 lt 5 lc 1 t "",  \
     'CO_rates_10AU.dat' u ($1/10):($2)      axis x1y1 w lines lw 3 lt 1 lc 2 t "10 AU",  \
     'CO_rates_10AU.dat' u ($1/10):($3/$2)   axis x1y2 w lines lw 3 lt 3 lc 2 t "",  \
     'CO_rates_10AU.dat' u ($1/10):($3)      axis x1y2 w lines lw 3 lt 5 lc 2 t "",  \
     'CO_rates_100AU.dat' u ($1/100):($2 )   axis x1y1 w lines lw 3 lt 1 lc 3 t "100 AU",  \
     'CO_rates_100AU.dat' u ($1/100):($3/$2) axis x1y2 w lines lw 3 lt 3 lc 3 t "",  \
     'CO_rates_100AU.dat' u ($1/100):($3)    axis x1y2 w lines lw 3 lt 5 lc 3 t ""
