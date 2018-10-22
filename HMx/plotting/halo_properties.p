reset

unset multiplot

#Font file
cmsy='/Users/Mead/Fonts/cmsy10.pfb'

#Terminal command
if(print==0){set term aqua; sun='sun'}
if(print==1){set term post enh col fontfile cmsy; set output 'halo_properties.eps'; sun='{/cmsy10 \014}'}

#Initial white space
print ''

#File to plot
file(z)=sprintf('data/properies_z%1.1f.dat',z)

#Mass axis
set log x
set xlabel 'M / h^{-1} M_{'.sun.'}'
set format x '10^{%T}'

#Align left and right margins
set lmargin 10
set rmargin 2

#Set the redshift
if(!exists('z')){z=0.0}
print 'Redshift z: ', z
print ''

set multiplot layout 2,2

#Radius
set log y
set ylabel 'R / h^{-1} Mpc'

set key top left

plot file(z) u 1:2 w l lw 2 ti 'Lagrangian radius',\
     file(z) u 1:3 w l lw 2 ti 'Virial radius'

set log y
set ylabel '{/Symbol n}'

plot file(z) u 1:4 w l lw 2 lc -1 noti

unset log y
set ylabel 'c'

plot file(z) u 1:5 w l lw 2 lc -1 noti

set log y
set ylabel '{/Symbol s}'

plot file(z) u 1:6 w l lw 2 lc -1 noti

unset multiplot

