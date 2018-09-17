reset

if(!exists("print")){print=0}
if(print==0){set term aqua dashed font ',14'}
if(print==1){set term post enh col sol; set output 'power.eps'}
if(print==2){set term pdfcairo; set output 'power.pdf'}

#This gnuplot script plots the output file 'power.dat' that is spat out of HMcode.
#Simply load up gnuplot (type gnuplot in the terminal) and then type "gnuplot>load 'plot.p'"
#The plot should then be the non-linear spectrum at 16 redshifts

# Initial white space
print ''

fac=4.*pi/(2.*pi)**3

kmin=1e-3
kmax=1e2

if(!exists("void")){void=0}

if(!exists("delta")){delta=1}
if(delta==0) {pmin=1e-3; pmax=1e5; plab='P(k) / (h^{-1} Mpc)^3'}
if(delta==1) {pmin=1e-8; pmax=1e4; plab='{/Symbol D}^2(k)'}

set log x
set xrange [kmin:kmax]
set xlabel 'k / (h Mpc^{-1})'
set mxtics 10

set log y
set yrange [pmin:pmax]
set ylabel plab
set format y '10^{%T}'

col='red'

unset colorbox

fpower='data/power.dat'
fvoid='data/power_1void.dat'

#Key stuff
if(delta==0) set key top right
if(delta==1) set key top left

# write
print 'delta = 0: Plot P(k)'
print 'delta = 1: Plot Delta^2(k)'
print 'delta = ', delta
print ''

print 'void = 0: Not plotting voids'
print 'void = 1: Plot voids'
print 'void = ', void
print ''

#Now do the actual plotting
if(delta==0){
if(void==0) {
plot fpower u 1:($2/(fac*$1**3)) w l lc -1       dt 1 lw 3 ti 'Linear',\
     fpower u 1:($3/(fac*$1**3)) w l lc rgb col  dt 2 lw 3 ti '2-halo',\
     fpower u 1:($4/(fac*$1**3)) w l lc rgb col  dt 3 lw 3 ti '1-halo',\
     fpower u 1:($5/(fac*$1**3)) w l lc rgb col  dt 1 lw 3 ti 'Total'
}

if(void==1) {
plot fpower u 1:($2/(fac*$1**3)) w l lc -1       dt 1 lw 3 ti 'Linear',\
     fpower u 1:($3/(fac*$1**3)) w l lc rgb col  dt 2 lw 3 ti '2-halo',\
     fpower u 1:($4/(fac*$1**3)) w l lc rgb col  dt 3 lw 3 ti '1-halo',\
     fvoid  u 1:($2/(fac*$1**3)) w l lc rgb col  dt 4 lw 3 ti '1-void',\
     fpower u 1:($5/(fac*$1**3)) w l lc rgb col  dt 1 lw 3 ti 'Total'
}
}

if(delta==1){
if(void==0) {
plot fpower u 1:2 w l lc -1       dt 1 lw 3 ti 'Linear',\
     fpower u 1:3 w l lc rgb col  dt 2 lw 3 ti '2-halo',\
     fpower u 1:4 w l lc rgb col  dt 3 lw 3 ti '1-halo',\
     fpower u 1:5 w l lc rgb col  dt 1 lw 3 ti 'Total'
}

if(void==1) {
plot fpower u 1:2 w l lc -1       dt 1 lw 3 ti 'Linear',\
     fpower u 1:3 w l lc rgb col  dt 2 lw 3 ti '2-halo',\
     fpower u 1:4 w l lc rgb col  dt 3 lw 3 ti '1-halo',\
     fvoid  u 1:2 w l lc rgb col  dt 4 lw 3 ti '1-void',\
     fpower u 1:5 w l lc rgb col  dt 1 lw 3 ti 'Total'
}
}




