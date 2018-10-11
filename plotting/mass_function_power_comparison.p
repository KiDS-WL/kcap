unset multiplot
reset session

ST='data/power_ShethTormen.dat'
Ti='data/power_Tinker.dat'

set log x
set format x ''

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set key top left

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

plot ST u 1:5 w l lc -1 lw 2 ti 'Sheth & Tormen (1999) mass function',\
     Ti u 1:5 w l lc 1  lw 2 ti 'Tinker et al. (2010) mass function'

set xlabel 'k / h Mpc^{-1}'
set format x

unset log y
set ylabel 'P_{Tinker}(k) / P_{Sheth-Tormen}(k)'
set yrange [*:1.02]
set format y

plot 1 w l lt -1 noti,\
     '<paste '.Ti.' '.ST.'' u 1:($5/$10) w l lw 2 lc 1 noti

unset multiplot
