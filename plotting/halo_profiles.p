reset
unset multiplot

cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(!exists('print')){print=0}
if(print==0){set term aqua dashed; msun='sun'}
if(print==1){set term post enh col font ',10' fontfile cmsy; msun='{/cmsy10 \014}'}

#File locations
file(m,z)=sprintf('diagnostics/halo_profile_m%i_z%1.1f.dat',m,z)
UPP(m,z)=sprintf('diagnostics/UPP/halo_profile_m%i_z%1.1f.dat',m,z)

#Initial white space
print ''

#Set the type of profiles to plot
if(!exists('type')){type='matter'}
print 'Options for *type*: matter, pressure'
print 'type: ', type
print ''

#Set the redshift
if(!exists('z')){z=0.0}
print 'Available z = 0.0, 0.5, 1.0, 2.0'
print 'Your z = ', z
print ''

m=13
print 'Example file: ', file(m,z)
print 'Example UPP: ', UPP(m,z)
print ''

#Title for individual plots
tits(m,z)=sprintf('M = 10^{%i} h^{-1} M_{'.msun.'}; z = %1.1f',m,z)

#Mass range for plot (10^m)
m1=13
m2=m1+2

if(type eq 'matter'){

if(print==1){set output 'halo_profiles.eps'}

if(!exists('pow')){pow=2}
if(pow==2){mtit=m2}
if(pow==3){mtit=m1}
#Power for y-axis 4*pi * r**pow * rho(r)
#pow==2 is good for linear x axis
#pow==3 is good for log x axis

#x axis
if(pow==2) {rmin=0.; rmax=1.1; unset log x}
if(pow==3) {rmin=1e-2; rmax=1.; set log x}
set xrange [rmin:rmax]
#set xrange [1e-2:1e1] #This is the range in the Fedeli paper
#set xlabel 'r / h^{-1} Mpc'
set xlabel 'r / r_v'

#y axis
rhomin=1e-3
if(pow==2){rhomax=10.}
if(pow==3){rhomax=1.}
set log y
set yrange [rhomin:rhomax]
set ylabel '4{/Symbol p} r^'.pow.' {/Symbol r}(r) / M    [h Mpc^{-1}]' offset 2
set mytics 10

set multiplot layout 1,3

do for [m=m1:m2] {

set title tits(m,z)

ti_CDM=''
ti_gas=''
ti_stars=''
ti_bound=''
ti_free=''
if(m==mtit) {ti_CDM='CDM'}
if(m==mtit) {ti_gas='Gas'}
if(m==mtit) {ti_stars='Stars'}
if(m==mtit) {ti_bound='Bound gas'}
if(m==mtit) {ti_free='Free gas'}

plot file(m,z) u 1:(4.*pi*($1**pow)*$2) w l lw 3 dt 1 lc rgb 'black' ti ti_CDM,\
     file(m,z) u 1:(4.*pi*($1**pow)*$3) w l lw 3 dt 1 lc rgb 'red'   ti ti_gas,\
     file(m,z) u 1:(4.*pi*($1**pow)*$4) w l lw 3 dt 1 lc rgb 'blue'  ti ti_stars,\
     file(m,z) u 1:(4.*pi*($1**pow)*$5) w l lw 3 dt 2 lc rgb 'red'   ti ti_bound,\
     file(m,z) u 1:(4.*pi*($1**pow)*$6) w l lw 3 dt 3 lc rgb 'red'   ti ti_free

}

unset multiplot

}

if(type eq 'pressure'){

if(print==1){set output 'halo_pressure.eps'}

if(!exists('ilog')){ilog=0}
if(ilog==0){rmin=0.; rmax=1.1; set xrange [rmin:rmax]}
if(ilog==1){rmin=1e-3; rmax=1e1; set log x; set xrange [rmin:rmax]}
#if(ih==0)   {set xlabel 'r / Mpc'}
set xlabel 'r / rv'

pmin=1e-4
pmax=1e1
if(ilog==0){set yrange [*:*]}#{set yrange [0:3.5]}
if(ilog==1){set log y; set yrange [pmin:pmax]; set format y '10^{%T}'; set mytics 10}
set ylabel '4{/Symbol p}r^2 P_e(r,M) / (eV cm^{-3} h^{-2} Mpc^2)'

set multiplot layout 1,3

do for [m=m1:m2] {

set title tits(m,z)

plot file(m,z) u 1:(4.*pi*$1*$1*$7) w l lw 3 dt 1 lc rgb 'red' ti 'Halo model',\
     UPP(m,z)  u 1:(4.*pi*$1*$1*$7) w l lw 3 dt 2 lc rgb 'black' ti 'UPP'


}

unset multiplot

}
