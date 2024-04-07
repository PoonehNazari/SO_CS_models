
set term postscript eps enhanced color 20 "Helvetica"
set out "C+_1AU.eps"
      
set logscale x
set xrange [0.1:1.000E+00]
set xtics 10
set mxtics 10
set xlabel "Radius (AU)"

set yrange [0:0.8] 
set ytics 0.1
set mytics 10
set ylabel "Z/R"

set logscale cb 
set cbrange [1e-12:0.0001]
set cbtics 10
set mcbtics 10
set format cb "10^{%T}"

set palette model RGB
set palette defined (0 "black", 1 "dark-violet", 2 "royalblue", 3 "sea-green", 4 "gold", 5 "dark-orange", 6 "dark-red")

set pm3d map
     
splot "C+_1AU.dat" u ($1):($2/$1):($3) t "" 

