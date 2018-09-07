reset

if(!exists('print')){print=0}
if(print==0) set term aqua
if(print==1) set term post enh col; set output 'kernels.eps'

k1='projection/kernel1.dat'
k2='projection/kernel2.dat'

set lmargin 15
set rmargin 15

set multiplot layout 2,1

set xlabel 'r / (h^{-1} Mpc)'

set ylabel 'X_1(r)'
set y2label 'X_2(r)'

set y2tics
set ytics nomirror
set format y2

plot k1 u 1:3 axes x1y1 w l lw 5 ti 'Kernel 1',\
     k2 u 1:3 axes x1y2 w l lw 3 ti 'Kernel 2'

#zmin=1e-2
#zmax=1e2
zmin=0.
zmax=10.
set xlabel 'z'
#set log x
set xrange [zmin:zmax]

plot k1 u 2:3 axes x1y1 w l lw 5 ti 'Kernel 1',\
     k2 u 2:3 axes x1y2 w l lw 3 ti 'Kernel 2'

unset multiplot
