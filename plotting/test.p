reset

set term aqua dashed

bench(name)=sprintf('benchmarks/power_HMcode_%s.txt',name)
power(name)=sprintf('data/power_HMx_%s_full.dat',name)

tests="'Mead' 'basic' 'standard'"

set ylabel 'P_{HMx}(k) / P_{HMcode}(k)'

set xlabel 'k / h^{-1} Mpc'
set log x

na=16

dx=1e-2
ddx=1e-3

set palette defined (1 'pink', 2 'grey')
unset colorbox

plot 1 w l lt -1 noti,\
     1+dx w l dt 2 lc -1 noti,\
     1-dx w l dt 2 lc -1 noti,\
     1+ddx w l dt 2 lc -1 noti,\
     1-ddx w l dt 2 lc -1 noti,\
     for [i=1:words(tests)] for [j=1:na] '<paste '.bench(word(tests,i)).' '.power(word(tests,i)).'' \
     u 1:(column(1+j)/column(2+j+na)):(real(j-1)/real(na-1)) w l lw 1 lc palette dt i noti
