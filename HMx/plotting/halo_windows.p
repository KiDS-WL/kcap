reset
unset multiplot

#Set terminal
cmsy='/Users/Mead/Fonts/cmsy10.pfb'
if(!exists('print')){print=0}
if(print==0){set term aqua dashed; sun='sun'}
if(print==1){set term post enh col fontfile cmsy; set output 'halo_windows.eps', sun='{/cmsy10 \014}'}

#file location
file(m,z)=sprintf('diagnostics/halo_window_m%i_z%1.1f.dat',m,z)

#x axis
kmin=1e-1
kmax=1e2
set xrange [kmin:kmax]
set xlabel 'kr_v'
set log x

#y axis
wmin=3.16e-3
wmax=3.16e0
set yrange [wmin:wmax]
set ylabel 'W_i(k,M) {/Symbol r}_m / M'
set log y

#Function for plot title
tits(m)=sprintf('M = 10^{%i} h^{-1} M_{'.sun.'}; z = 0',m)

#Fix the redshift
z=0.0

#Set the mass range in log(M)
m1=13
m2=m1+2

#Set multiplot mode
set multiplot layout 1,3

#Loop over masses to make three plots
do for [m=m1:m2] {

#Get the individial plot titles
set title tits(m)

#Make the individual plot
plot file(m,z) u 1:2 w l dt 1 lw 3 lc rgb 'black' ti 'CDM',\
     file(m,z) u 1:3 w l dt 1 lw 3 lc rgb 'red' ti 'Gas',\
     file(m,z) u 1:4 w l dt 1 lw 3 lc rgb 'blue' ti 'Stars',\
     file(m,z) u 1:5 w l dt 2 lw 3 lc rgb 'red' ti 'Bound gas',\
     file(m,z) u 1:6 w l dt 3 lw 3 lc rgb 'red' ti 'Free gas'

}

#Unset multiplot
unset multiplot
