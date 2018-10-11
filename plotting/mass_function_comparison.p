unset multiplot
reset

if(!exists('print')){print=0}
if(print==0){set term aqua; sun='sun'}

print ''

bnu='data/bnu_functions.dat'
gnu='data/gnu_functions.dat'
nm='data/mass_functions.dat'

if(!exists('iplot')){iplot=1}
print 'iplot = 1: g(nu) and b(nu)'
print 'iplot = 2: n(M)'
print 'iplot = 3: g(M) and b(M)'
print 'iplot = ', iplot
print ''

if(iplot==1 || iplot==3){

set log y

set lmargin 10
set rmargin 2

set multiplot layout 2,1

do for [i=1:2]{

if(iplot==1){

unset log x
set xlabel '{/Symbol n}'
set format x
set xrange [*:*]

if(i==1){file=gnu; set ylabel 'g({/Symbol n})'; set format y '10^{%T}'}
if(i==2){file=bnu; set key top left; set ylabel 'b({/Symbol n})'; set format y}

plot file u 1:3 w l lw 3 ti 'Press & Schecter (1974)',\
     file u 1:4 w l lw 3 ti 'Sheth & Tormen (1999)',\
     file u 1:5 w l lw 3 ti 'Tinker et al. (2010)'

}

if(iplot==3){

print 'Not sure that plotting g(nu) as g(M) with M(nu) is a reasonable thing to do'

mmin=1e8
mmax=1e16
set log x
set xlabel 'M / h^{-1} M_{'.sun.'}'
set format x '10^{%T}'
set xrange [mmin:mmax]

if(i==1){file=gnu; set ylabel 'g(M)'; set format y '10^{%T}'}
if(i==2){file=bnu; set key top left; set ylabel 'b(M)'; set format y}

plot file u 2:3 w l lw 3 ti 'Press & Schecter (1974)',\
     file u 2:4 w l lw 3 ti 'Sheth & Tormen (1999)',\
     file u 2:5 w l lw 3 ti 'Tinker et al. (2010)'

}

}

}

if(iplot==2){

mmin=1e8
mmax=1e16
set xlabel 'M / h^{-1} M_{'.sun.'}'
set log x
set xrange [mmin:mmax]

set ylabel 'M^2 n(M) / {/Symbol r}_m'
set log y
set format y '10^{%T}'

plot nm u 2:5 w l lw 3 noti

}

print ''

unset multiplot


