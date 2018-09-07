reset

if(print==0) set term aqua dashed
if(print==1) set term post enh col sol; set output 'linear_ratio.eps'

#Initial white space
print ''

#k axis
kmin=1e-3
kmax=1e0
set log x
set xrange [kmin:kmax]
set xlabel 'k / (h Mpc^{-1})'
set mxtics 10

#Ratio axis
dp=0.1
set yrange [1-dp:1+dp]
set ylabel 'P(k) / P_{lin}(k)'

col='red'

#File to plot
fpower='data/power.dat'
print 'Plotting: ', fpower
print ''

#Key stuff
set key bottom left

#Now do the actual plotting
plot 1.00 w l lc -1 dt 1 noti,\
     1.01 w l lc -1 dt 2 noti,\
     0.99 w l lc -1 dt 2 noti,\
     1.05 w l lc -1 dt 3 noti,\
     0.95 w l lc -1 dt 3 noti,\
     fpower u 1:($3/$2) w l lc rgb col  dt 2 lw 3 ti '2-halo',\
     fpower u 1:($5/$2) w l lc rgb col  dt 1 lw 3 ti 'Total'




