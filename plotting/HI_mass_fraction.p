reset

sun='sun'

set log x
set xlabel 'M / h^{-1} M_{'.sun.'}'
set format x '10^{%T}'

set ylabel 'M_{HI} / M'

nz=6

zs="'z = 0' 'z = 1' 'z = 2' 'z = 3' 'z = 4' 'z = 5'"

plot for [i=1:nz] 'data/HI_mass_fraction.dat' u 1:(column(i+1)) w l lw 3 lc i ti word(zs,i)
