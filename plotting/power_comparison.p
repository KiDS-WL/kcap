unset multiplot
reset

#This gnuplot script plots the output file 'power.dat' that is spat out of HMcode.
#Simply load up gnuplot (type gnuplot in the terminal) and then type "gnuplot>load 'plot.p'"
#The plot should then be the non-linear spectrum at 16 redshifts

if(!exists('print')){print=0}
if(print==0){set term aqua dashed}
if(print==1){set term post enh col sol; set output 'power_comparison.eps'}

kmin=1e-3
kmax=1e2
set log x
set xrange [kmin:kmax]

pmin=1e-7
pmax=1e4
set log y
set yrange [pmin:pmax]
#set yrange [*:*]
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

#File to compare against
if(!exists('com')){com='benchmarks/power_HMcode_Mead.txt'}

new='data/power_hm.dat'
p1h='data/power_1h.dat'
p2h='data/power_2h.dat'

print ''
print 'Comparing power_hm.dat to: ', com
print 'Set variable *com* as comparison file'
print ''

set key top left

col_red(i)=sprintf("#%1x%1x0000",i-1,i-1)
col_green(i)=sprintf("#00%1x%1x00",i-1,i-1)

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

nz=16

plot for[i=1:nz]  com u 1:(column(i+1)) w l lw 2 dt 1 lc rgb col_red(i) noti,\
     for[i=1:nz]  new u 1:(column(i+1)) w l lw 2 dt 1 lc rgb col_green(i) noti,\
     for[i=nz:nz] p2h u 1:(column(i+1)) w l lw 2 dt 2 lc rgb col_green(i) noti,\
     for[i=nz:nz] p1h u 1:(column(i+1)) w l lw 2 dt 3 lc rgb col_green(i) noti

set xlabel 'k / (h Mpc^{-1})'
set format x

unset log y
set yrange [*:*]
#set yrange [0.98:1.02]
#set yrange [0:2]
set format y
set ylabel 'P_{new}(k) / P_{old}(k)'

dy=0.01
ddy=0.001

plot 1 w l ls -1 noti,\
     1.-ddy w l lw 1 lc -1 dt 3 ti '0.1%',\
     1.+ddy w l lw 1 lc -1 dt 3 noti,\
     1.-dy  w l lw 1 lc -1 dt 2 ti '1%',\
     1.+dy  w l lw 1 lc -1 dt 2 noti,\
     for[i=1:nz] '<paste '.new.' '.com.'' u 1:(column(i+1)/column(i+1+17)) w l lw 2 dt 1 lc rgb col_green(i) noti
     #for[i=1:nz] '<paste '.new.' '.com.'' u 1:(column(i+1)/column(35-i)) w l lw 2 dt 1 lc rgb col_green(i) noti

unset multiplot



