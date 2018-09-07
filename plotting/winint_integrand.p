reset

file='winint/integrand.dat'

set xlabel 'r / r_v'

set ylabel '4{/Symbol p}r^2 sinc(kr) {/Symbol r}(r)'

kmin=1e-1
kmax=1e2
nk=16
set cbrange [kmin:kmax]
set cblabel 'k / (h Mpc^{-1})'
set palette defined (1 'black', 2 'blue', 3 'light-blue')
set log cb

plot 0 w l lt -1 noti,\
     for [i=1:nk] file u 1:(column(i+1)):(exp(log(kmin)+log(kmax/kmin)*real(i-1)/real(nk-1))) w l lw 3 lc palette noti
