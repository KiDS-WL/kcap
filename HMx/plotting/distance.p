reset

#File location
file='projection/distance.dat'

#x axis
zmin=0.
zmax=4.
set xlabel 'z'
set xrange [zmin:zmax]

#y axis
rmin=0.
rmax=6000.
set ylabel 'r / (h^{-1} Mpc)'
set yrange [rmin:*]

#Set the key
set key top left

#Make the plot
plot file u 1:2 w p pt 7 ps 1 noti,\
     file u 1:2 w l lw 3 ti 'Comoving distance',\
     file u 1:3 w l lw 3 ti 'Comoving angular distance'

