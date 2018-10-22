reset

sun='sun'

labels="'L = infinity' 'L = 2048 Mpc/h' 'L = 1024 Mpc/h' 'L = 512 Mpc/h' 'L = 256 Mpc/h' 'L = 128 Mpc/h' 'L = 64 Mpc/h' 'L = 32 Mpc/h'"

nbox=8

set xlabel '{/Symbol n}'

set ylabel 'M / h^{-1} M_{'.sun.'}'
set log y
set yrange [1e10:1e16]

set key bottom right

plot for [i=1:nbox] 'data/nu_mass_Lbox.dat' u 1:(column(i+1)) w l lw 3 ti word(labels,i)
