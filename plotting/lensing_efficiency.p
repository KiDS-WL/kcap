unset multiplot
reset

if(print==0) set term aqua dashed
if(print==1) set term post enh col; set output 'lensing_efficiency.eps'

file='data/lensing_efficiency.dat'

set ylabel 'q(r)'
set yrange [0:1]

set multiplot layout 2,1

set xlabel 'r / (Mpc/h)'

plot file u 1:3 w l lw 3 noti

set xlabel 'z'

plot file u 2:3 w l lw 3 noti

unset multiplot
