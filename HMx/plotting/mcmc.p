reset

np=12

plot for [i=1:np] 'data/mcmc.dat' u (column(0)):(column(i+1)) w p noti
