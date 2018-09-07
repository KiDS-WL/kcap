unset multiplot
reset

if(!exists('print')){print=0}
if(print==0){set term aqua dashed}

# Initial white space
print ''

#Alonso(i,z)=sprintf('/Users/Mead/Physics/CCL/tests/benchmark/pk_hm_c%d_z%d.txt',i,z)
Alonso(i,z)=sprintf('/Users/Mead/Physics/CCL_tests/pk_hm_c%d_z%d.txt',i,z)
#Alonso(i,z)=sprintf('benchmarks/Alonso_set1/pk_hm_c%d_z%d.txt',i,z)
#Alonso(i,z)=sprintf('benchmarks/Alonso_set2/pk_hm_c%d_z%d.txt',i,z)
HMx(i,z)=sprintf('benchmarks/HMx_power_model%d_z%d.txt',i,z)

#kmin=1e-3
#kmax=1e1
set log x
#set xrange [kmin:kmax]
set xrange [*:*]

set log y
set format y '10^{%T}'
set ylabel '{/Symbol D}^2(k)'

f=(4.*pi)/(2.*pi)**3

set key top left

if(!exists('itype')){itype=1}
print 'itype = 1: Full power'
print 'itype = 2: Linear theory'
print 'itype = 3: Two-halo'
print 'itype = 4: One-halo'
#if(itype==1){c1=5; c2=2; print 'itype ', itype, ': Full power'}
#if(itype==2){c1=2; c2=5; print 'itype ', itype, ': Linear power'}
#if(itype==3){c1=3; c2=4; print 'itype ', itype, ': Two-halo term'}
#if(itype==4){c1=4; c2=3; print 'itype ', itype, ': One-halo term'}
if(itype==1){c1=5; c2=5; print 'itype ', itype, ': Full power'}
if(itype==2){c1=2; c2=2; print 'itype ', itype, ': Linear power'}
if(itype==3){c1=3; c2=3; print 'itype ', itype, ': Two-halo term'}
if(itype==4){c1=4; c2=4; print 'itype ', itype, ': One-halo term'}
print ''

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

z=0.

plot for [i=1:3] HMx(i,z) u 1:(column(c1)/10**i) w l lw 2 lc i dt 2 ti 'HMx',\
     for [i=1:3] Alonso(i,z) u 1:((f*$1**3)*column(c2)/10**i) w l lw 2 lc i dt 1 ti 'Alonso'

set xlabel 'k / h Mpc^{-1}'
set format x

unset log y
#set yrange [0.99:1.01]
set format y

plot 1 w l lc -1 dt 1 noti,\
     1.01 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.001 w l lc -1 dt 2 noti,\
     0.999 w l lc -1 dt 2 noti,\
     for [i=1:3] '<paste '.HMx(i,z).' '.Alonso(i,z).'' u 1:(column(c1)/((f*$1**3)*column(5+c2))) w l lw 2 lc i dt 1 noti

unset multiplot
