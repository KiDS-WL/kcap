reset

pow_norm='data/power_hm.dat'
pow_scat='data/power_scatter_hm.dat'

set log x
set xlabel 'k / h^{-1} Mpc'

set ylabel 'P_{scatter}(k) / P_{normal}(k)'

na=16
amin=0.1
amax=1.
set cbrange [amin:amax]
set palette cubehelix negative
set cblabel 'a'

plot 1 w l lt -1 noti,\
     for [i=1:na] '<paste '.pow_scat.' '.pow_norm.'' u 1:(column(1+i)/column(2+na+i)):(amin+(amax-amin)*(real(i-1)/real(na-1))) w l lc palette lw 3 noti
