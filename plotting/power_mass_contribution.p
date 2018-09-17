reset

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(!exists('print')){print=0}
if(print==0){set term aqua dashed size 1000,500; sun='sun'}
if(print==1){set term post enh col fontfile cmsy ',10' size 10,5; set output 'paper/power_mass_contribution.eps'; sun='{/cmsy10 \014}'}

# data file
power(f1,f2,m)=sprintf('data/power_%d%d_m%d.dat',f1,f2,m)

# Mass range
m1=10
m2=16

# wavenumber axis
kmin=1e-2
kmax=1e1
set log x
set xrange [kmin:kmax]
set xlabel 'k / h Mpc^{-1}'

set multiplot layout 1,2

# power axis
dmin=1e-3
dmax=1e3
set log y
set yrange [dmin:dmax]
set ylabel '{/Symbol D}^2_{mm}(k)'
set format y '10^{%T}'

set palette defined ( 0 "light-blue", 1 "blue", 2 "black" )
set cblabel 'log_{10} (M / h^{-1} M_{'.sun.'})'

plot for [i=m1:m2] power(0,0,i) u 1:5:(i) w l lw 3 dt 1 lc palette noti#,\
     for [i=m1:m2] power(0,0,i) u 1:3:(i) w l lw 3 dt 2 lc palette noti,\
     for [i=m1:m2] power(0,0,i) u 1:4:(i) w l lw 3 dt 3 lc palette noti

# power axis
dmin=1e-7
dmax=1e-1
set log y
set yrange [dmin:dmax]
set ylabel '{/Symbol D}^2_{mp}(k) / eV cm^{-3}'
set format y '10^{%T}'

set palette defined ( 0 "pink", 1 "red", 2 "black" )
set cblabel 'log_{10} (M / h^{-1} M_{'.sun.'})'

plot for [i=m1:m2] power(0,6,i) u 1:5:(i) w l lw 3 dt 1 lc palette noti#,\
     for [i=m1:m2] power(0,6,i) u 1:3:(i) w l lw 3 dt 2 lc palette noti,\
     for [i=m1:m2] power(0,6,i) u 1:4:(i) w l lw 3 dt 3 lc palette noti

unset multiplot
