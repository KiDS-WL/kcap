reset
unset multiplot

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(print==0){set term aqua; ell='l'}
if(print==1){set term post enh col fontfile cmmi; ell='{/cmmi10 \140}'; set output 'xcoef.eps'}

Cl_11='data/cl_first.dat'
Cl_22='data/cl_second.dat'
Cl_12='data/cl_hm.dat'

ellmin=1e1
ellmax=1e4
set log x
set xrange [ellmin:ellmax]

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

set log y
set ylabel ''.ell.'('.ell.'+1)C('.ell.') / 2{/Symbol p}'
set format y '10^{%T}'

set key bottom right

plot Cl_11 u 1:3 w l lw 3 ti 'Autospectra: field 1',\
     Cl_22 u 1:3 w l lw 3 ti 'Autospectra: field 2',\
     Cl_12 u 1:3 w l lw 3 ti 'Cross spectra'

set xlabel ''.ell.''
set format x '10^{%T}'

rmin=0.
rmax=1.05
unset log y
set ylabel 'r('.ell.')'
set yrange [rmin:rmax]
set format y

plot 1 w l lt -1 noti,\
     '<paste '.Cl_12.' '.Cl_11.' '.Cl_22.'' u 1:($2/sqrt($5*$8)) w l lc -1 lw 3 noti

unset multiplot
