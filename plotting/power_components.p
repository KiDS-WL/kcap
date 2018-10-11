reset

# Terminal options
if(!exists('print')){print=0}
if(print==0){set term aqua dashed}
if(print==1){set term post enh col font ',12' size 7,7; set output 'paper/power_components.eps'}

# Initial white space
print ''

# File name
power(z,f1,f2)=sprintf('data/power_z%1.1f_%i%i.dat',z,f1,f2)

# k axis
kmin=1e-2
kmax=1e1
set xrange [kmin:kmax]
set log x
set xlabel 'k / h Mpc^{-1}'

# Delta^2 axis
dmin=1e-10
dmax=1e3
set yrange [dmin:dmax]
set log y
set format y '10^{%T}'
set ylabel '{/Symbol D}^2_{uv}(k)'

# Pressure field multiplications
f1=sqrt(10)
f2=f1**2
print 'Pressure field multiplied by: ', f1
print 'Pressure auto-spectrum multiplied by: ', f2

# Set the redshift
if(!exists('z')){z=0.}
print 'z = ', z

# Final white space
print ''

# Key position
set key top left

set label sprintf('z = %1.1f', z) at screen 0.9,0.15

# Plot
plot NaN w l lw 3 lc -1 dt 2 ti 'Two-halo term',\
     NaN w l lw 3 lc -1 dt 3 ti 'One-halo term',\
     NaN w l lw 3 lc -1 dt 1 ti 'Halo model power',\
     power(z,0,0) u 1:5       w l lc 1 dt 1 lw 3 ti 'Matter',\
     power(z,1,1) u 1:3       w l lc 2 dt 2 lw 3 noti,\
     power(z,1,1) u 1:4       w l lc 2 dt 3 lw 3 noti,\
     power(z,1,1) u 1:5       w l lc 2 dt 1 lw 3 ti 'CDM',\
     power(z,2,2) u 1:3       w l lc 3 dt 2 lw 3 noti,\
     power(z,2,2) u 1:4       w l lc 3 dt 3 lw 3 noti,\
     power(z,2,2) u 1:5       w l lc 3 dt 1 lw 3 ti 'Gas',\
     power(z,3,3) u 1:3       w l lc 4 dt 2 lw 3 noti,\
     power(z,3,3) u 1:4       w l lc 4 dt 3 lw 3 noti,\
     power(z,3,3) u 1:5       w l lc 4 dt 1 lw 3 ti 'Stars',\
     power(z,6,6) u 1:(f2*$3) w l lc 6 dt 2 lw 3 noti,\
     power(z,6,6) u 1:(f2*$4) w l lc 6 dt 3 lw 3 noti,\
     power(z,6,6) u 1:(f2*$5) w l lc 6 dt 1 lw 3 ti 'Electron pressure'

unset label
