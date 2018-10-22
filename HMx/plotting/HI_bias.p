unset multiplot
reset

if(!exists('print')){print=0}
if(print==0){set term aqua}
if(print==1){set term post enh col; set output 'HI_bias.eps'}

mm='data/power_mm_hm.dat'
mf='data/power_mf_hm.dat'
ff='data/power_ff_hm.dat'
sim(f1,f2,z)=sprintf('/Users/Mead/Physics/people/Villaescusa_Navarro/HI_bias/b_%s%s_mm_z=%d.txt',f1,f2,z)

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

set multiplot layout 1,2

set title 'Bias from auto-correlation'

plot 1 w l lt -1 noti,\
     for [z=0:5]  sim('HI','HI',z) w p pt 7 ps .5 lc z+1 noti,\
     for [i=1:nz] '<paste '.ff.' '.mm.'' u 1:(sqrt(column(1+i)/column(nz+2+i))) w l lc i dt 1 lw 3 ti word(zs,i)


set title 'Bias from cross-correlation'

plot 1 w l lt -1 noti,\
     for [z=0:5]  sim('HI','m',z) w p pt 7 ps .5 lc z+1 noti,\
     for [i=1:nz] '<paste '.mf.' '.mm.'' u 1:(column(1+i)/column(nz+2+i)) w l lc i dt 1 lw 3 noti

unset multiplot
