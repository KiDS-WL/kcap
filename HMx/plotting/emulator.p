unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col font ',12'; set output 'emulator.eps'}

# Initial white space
print ''

file(n,z)=sprintf('data/cosmo%d_z%d.dat',n,z)

# Number of redshifts
nz=4

if(!exists('n')){n=10}
print 'Plotting *n* cosmolgies:', n

kmin=1e-2
kmax=1e1
set log x
set xlabel ''
set format x ''
set xrange [kmin:kmax]

dy=0.12
ddy=0.05
set yrange [1.-dy:1+dy]
set ylabel
set ylabel 'P(k) / P_{emu}(k)'

set lmargin 10
set rmargin 2

y2=0.98
y1=0.1
dy=(y2-y1)/4

# Final white space
print ''

set multiplot layout 4,1

do for [j=1:nz]{

if(j==1){set tmargin at screen y2;      set bmargin at screen y2-dy; set label 'z = 0.0' at graph 0.93,0.9}
if(j==2){set tmargin at screen y2-dy;   set bmargin at screen y2-2*dy; set label 'z = 0.5' at graph 0.93,0.9}
if(j==3){set tmargin at screen y1+2*dy; set bmargin at screen y1+dy; set label 'z = 1.0' at graph 0.93,0.9}
if(j==4){set tmargin at screen y1+dy;   set bmargin at screen y1; set label 'z = 2.0' at graph 0.93,0.9}
if(j==4){set xlabel 'k / h Mpc^{-1}'; set format x}

plot 1 w l lt -1 noti,\
     1.+ddy w l lc -1 dt 2 noti,\
     1.-ddy w l lc -1 dt 2 noti,\
     for [i=0:n] file(i,j) u 1:($2/$3) w l noti

unset label

}

unset multiplot
