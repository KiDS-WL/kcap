reset

if(!exists('print')){print=0}
if(print==0){set term aqua}
if(print==1){set term post enh col; set output 'HI_bias.eps'}

mm='data/power_mm_full.dat'
mf='data/power_mf_full.dat'
ff='data/power_ff_full.dat'

kmin=0.1
kmax=70.
set xrange [kmin:kmax]
set log x
set xlabel 'k / h Mpc^{-1}'

bmin=0.3
bmax=20.
set ylabel 'b_{HI}(k)'
set yrange [bmin:bmax]
set log y

# Total number of redshifts
nz=6

zs="'z=0' 'z=1' 'z=2' 'z=3' 'z=4' 'z=5'"

set key top left

plot 1 w l lt -1 noti,\
     NaN w l lw 3 lc -1 dt 1 ti 'Bias from auto-correlation',\
     NaN w l lw 3 lc -1 dt 2 ti 'Bias from cross-correlation',\
     for [i=1:nz] '<paste '.ff.' '.mm.'' u 1:(sqrt(column(1+i)/column(nz+2+i))) w l lc i dt 1 lw 3 ti word(zs,i),\
     for [i=1:nz] '<paste '.mf.' '.mm.'' u 1:(column(1+i)/column(nz+2+i))       w l lc i dt 2 lw 3 noti
