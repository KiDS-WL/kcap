unset multiplot
reset

set term aqua dashed size 1000,800 font ',16'

bench(name)=sprintf('benchmarks/power_HMcode_%s.txt',name)
power(name)=sprintf('data/power_HMx_%s_hm.dat',name)

print ''
print 'itest = 1: Test against HMcode'
print 'itest = 2: Test hydro'
if(!exists('itest')){itest=1}
print 'itest = ', itest
print ''

set xlabel 'k / h^{-1} Mpc'
set log x

dy=1e-2
ddy=1e-3

set palette defined (1 'red', 2 'pink', 3 'grey')
unset colorbox

if(itest==1){

na=16

set ylabel 'P_{HMx}(k) / P_{HMcode}(k)'

tests="'Mead' 'basic' 'standard'"

plot 1 w l lt -1 noti,\
     1+dy w l dt 2 lc -1 noti,\
     1-dy w l dt 2 lc -1 noti,\
     1+ddy w l dt 2 lc -1 noti,\
     1-ddy w l dt 2 lc -1 noti,\
     for [i=1:words(tests)] for [j=1:na] '<paste '.bench(word(tests,i)).' '.power(word(tests,i)).'' \
     u 1:(column(1+j)/column(2+j+na)):(real(j-1)/real(na-1)) w l lw 2 lc palette dt i noti

}

if(itest==2){

bench(z,i,j)=sprintf('benchmarks/power_z%1.1f_%d%d.txt',z,i,j)
test(z,i,j)=sprintf('data/power_z%1.1f_%d%d.dat',z,i,j)

set ylabel 'P_{HMx}(k) / P_{benchmark}(k)'

fields="'matter' 'CDM' 'gas' 'stars' 'electron pressure'"
array f[5]
f[1]=0
f[2]=1
f[3]=2
f[4]=3
f[5]=6

array c[5]
c[1]=1
c[2]=2
c[3]=3
c[4]=4
c[5]=6

set multiplot layout 2,2

do for [iz=1:4]{

if(iz==1){z=0.0}
if(iz==2){z=0.5}
if(iz==3){z=1.0}
if(iz==4){z=2.0}

set label sprintf('z = %1.1f',z) at graph 0.1,0.9

plot 1 w l lt -1 noti,\
     1+dy w l dt 2 lc -1 noti,\
     1-dy w l dt 2 lc -1 noti,\
     1+ddy w l dt 2 lc -1 noti,\
     1-ddy w l dt 2 lc -1 noti,\
     for [i=1:5] for [j=i:5] '<paste '.test(z,f[i],f[j]).' '.bench(z,f[i],f[j]).'' u 1:($5/$10) w l lw 2 lc c[i] dt j noti

unset label

}

unset multiplot

}
