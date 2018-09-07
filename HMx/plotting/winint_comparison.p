unset multiplot
reset

set term aqua dashed

#kmin=1e-3
#kmax=1e3
set log x
#set xrange [kmin:kmax]
set mxtics 10

file(i)=sprintf('winint/results_%d.dat',i)

set key top left

set lmargin 10
set rmargin 2

n=7
do for [j=1:n] {

set multiplot layout 2,1

set xlabel ''
set format x ''

wmin=1e-5
wmax=2.
set log y
set ylabel 'W(k)'
set yrange [wmin:wmax]
set format y '10^{%T}'

set key bottom left

plot for [i=1:j] file(i) u 1:(+$2) w l dt 1 lc i lw 3 ti file(i),\
     for [i=1:j] file(i) u 1:(-$2) w l dt 2 lc i lw 3 noti

set format x
set xlabel 'k / (h Mpc^{-1})'
set format x '10^{%T}'

unset log y
set ylabel 'W_{i}(k) / W_{base}(k)'
set yrange [0.99:1.01]
set format y

error=1e-3
plot 1.+error w l dt 2 lc -1 noti,\
     1.-error w l dt 2 lc -1 noti,\
     for [i=1:j] '<paste '.file(i).' '.file(1).'' u 1:($2/$4) w l dt 1 lc i lw 3 noti

unset multiplot

if(j<n) {pause 0.5}

}
