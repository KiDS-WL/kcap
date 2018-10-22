unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed}
if(print==2) {set term qt}

# Initial white space
print ''

# Data files
data(base,type,icosmo,if1,if2,iz)=sprintf('%s_%s_cos%d_%d%d_z%d.dat',base,type,icosmo,if1,if2,iz)

# Set number of cosmoloies
if(!exists('ncos')){ncos=1}
print 'Number of cosmologies: ncos: ', ncos

if(!exists('nf')){nf=1}
print 'Number of fields: nf: ', nf

# Set number of redshifts
if(!exists('nz')){nz=1}
print 'Number of redshifts: nz: ', nz

if(!exists('base')){base='fitting/test'}
print 'Data file base: base: ', base
print 'Example data file: ', data(base,'best',1,1,1,1)
print ''

kmin=0.1
kmax=10.
set xrange [kmin:kmax]
set log x

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set key top left

set xlabel ''
set format x ''

set log y
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'

col(i,ni,j,nj,k,nk,l,nl)=nj*nk*nl*(i-1)+nk*nl*(j-1)+nl*(k-1)+l

plot for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:3 w p lc col(i,ncos,j,nz,j1,nf,j2,nf) dt 1 lw 2 noti 'BAHAMAS',\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:2 w l lc 0                            dt j lw 2 noti 'HMx',\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:2 w l lc col(i,ncos,j,nz,j1,nf,j2,nf) dt j lw 2 noti 'HMx'

set xlabel 'k / h^{-1} Mpc'
set format x

rmin=0.8
rmax=1.2
unset log y
set yrange [rmin:rmax]
set ylabel 'P_{HM}(k) / P_{sim}(k)'
set format y

plot 1 w l lt -1 noti,\
     0.95 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.01 w l lc -1 dt 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'orig',i,j1,j2,j) u 1:($2/$3) w l lc 0                            dt j lw 2 noti,\
     for [i=1:ncos] for [j=1:nz] for [j1=1:nf] for [j2=j1:nf] data(base,'best',i,j1,j2,j) u 1:($2/$3) w l lc col(i,ncos,j,nz,j1,nf,j2,nf) dt j lw 2 noti

unset multiplot
