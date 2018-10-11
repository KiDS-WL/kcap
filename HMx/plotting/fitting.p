unset multiplot
reset

set term aqua dashed

# Initial white space
print ''

# Data files
data(name,icosmo,iz)=sprintf('fitting/%scosmo%d_z%d.dat',name,icosmo,iz)

# Set number of cosmoloies
if(!exists('ncos')){ncos=1}
print 'Number of cosmologies: *ncos*: ', ncos

# Set number of redshifts
if(!exists('nz')){nz=1}
print 'Number of redshifts: *nz*: ', nz

set log x

set lmargin 10
set rmargin 2

# Final white space
print ''

set multiplot layout 2,1

set key top left

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

plot for [i=1:ncos] for [j=1:nz] data('',i,j)          u 1:3 w p lc i pt 1 noti 'BAHAMAS',\
     for [i=1:ncos] for [j=1:nz] data('original_',i,j) u 1:2 w l lc 0 dt j lw 2 noti 'HMx',\
     for [i=1:ncos] for [j=1:nz] data('',i,j)          u 1:2 w l lc i dt j lw 2 noti 'HMx'

set xlabel 'k / h^{-1} Mpc'
set format x

unset log y
set ylabel 'P_{HM}(k) / P_{sim}(k)'
set format y

plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.01 w l lc -1 dt 2 noti,\
     for [i=1:ncos] for [j=1:nz] data('original_',i,j) u 1:($2/$3) w l lc 0 dt j lw 2 noti,\
     for [i=1:ncos] for [j=1:nz] data('',i,j)          u 1:($2/$3) w l lc i dt j lw 2 noti

unset multiplot
