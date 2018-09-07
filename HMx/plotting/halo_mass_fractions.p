reset

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(!exists('print')){print=0}
if(print==0){set term aqua dashed; sun='sun'}
if(print==1){set term post enh col fontfile cmsy; set output 'halo_mass_fractions.eps'; sun='{/cmsy10 \014}'}

file='diagnostics/mass_fractions.dat'

set size square

set log x
set xlabel 'M / h^{-1} M_{'.sun.'}'
set xrange [1e10:1e16]
set format x '10^{%T}'
set mxtics 10

set log y
set ylabel 'Mass fraction'
#set yrange [3.162e-3:3.162e0]
set yrange [3.162e-3:1e0]
set format y

set key outside left box

om_b=0.05
om_m=0.3

plot om_b/om_m  w l lw 3 dt 2 lc rgb 'black' ti 'Universal baryon',\
     file u 1:2 w l lw 3 dt 1 lc rgb 'black' ti 'CDM',\
     file u 1:3 w l lw 3 dt 1 lc rgb 'red' ti 'Gas',\
     file u 1:4 w l lw 3 dt 1 lc rgb 'blue' ti 'Stars',\
     file u 1:5 w l lw 3 dt 2 lc rgb 'red' ti 'Bound gas',\
     file u 1:6 w l lw 3 dt 3 lc rgb 'red' ti 'Free gas'

