unset multiplot
reset

if(!exists('print')){print=0}
if(print==0){set term aqua}
if(print==1){set term post enh col; set output 'cored.eps'}

rcore_min=0.
rcore_max=0.1
ncore=16

power(i)=sprintf('data/power_cored_%d.dat',i)

set log x
set xrange [*:*]

set key top left

top=0.97
bot=0.1
dy=0.02
lef=0.1
rig=0.85

set lmargin at screen lef
set rmargin at screen rig

set multiplot layout 2,1

set tmargin at screen top
set bmargin at screen (top+bot+dy)/2.

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

#set palette cubehelix
set cblabel 'r_{core} / h^{-1} Mpc'
set colorbox vertical user origin rig+0.02, bot size 0.04,top-bot

plot for [i=1:ncore] power(i) u 1:5:(rcore_min+(rcore_max-rcore_min)*(i-1)/(ncore-1)) w l lw 3 lc palette noti

set tmargin at screen (top+bot-dy)/2.
set bmargin at screen bot

set xlabel 'k / h^{-1} Mpc'
set format x

rmin=0.
rmax=1.1
unset log y
set ylabel 'P(k) / P_{NFW}(k)'
set yrange [rmin:rmax]
set format y

unset colorbox

unset key

plot 1 w l lt -1 noti,\
     for [i=1:ncore] '<paste '.power(i).' '.power(1).'' u 1:($5/$10):(rcore_min+(rcore_max-rcore_min)*(i-1)/(ncore-1)) w l lw 3 lc palette noti

unset multiplot
