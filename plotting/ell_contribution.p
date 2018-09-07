unset multiplot
reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(print==0){set term aqua font ',10'; ell='l'}
if(print==1){set term post enh col fontfile cmmi ',10'; set output 'ell_contribution.eps'; ell='{/cmmi10 \140}'}

file(n)=sprintf('Limber/Cell_contrib_ell_%d.dat', n)

#Fix the y axis to have no label
set ylabel ''
set format y ''
set yrange [0.:*]

#Do the colour bar
n=16
ellmin=2**0
ellmax=2**(n-1)
set log cb
set cbrange [ellmin:ellmax]
set cblabel ''.ell.''
unset colorbox

#Get the plot dimensions
lef=0.02
rig=0.9
top=0.98
bot=0.1
ny=3
gap=0.07
dy=(top-bot-((ny-1)*gap))/ny

#Set the left and right margins
set lmargin at screen lef
set rmargin at screen rig

set multiplot

#Function of k plot
set tmargin at screen top
set bmargin at screen top-dy
kmin=1e-3
kmax=1e2
set xlabel 'k / (h Mpc^{-1})'
set log x
set xrange [kmin:kmax]
plot for [i=1:n] file(i) u 1:2:(2**(i-1)) w l lw 3 lc palette noti

#Function of z plot
set tmargin at screen top-dy-gap
set bmargin at screen top-2*dy-gap
zmin=1e-3
zmax=1e1
set xlabel 'z'
set log x
set xrange [zmin:zmax]
plot for [i=1:n] file(i) u 3:4:(2**(i-1)) w l lw 3 lc palette noti

#Function of R plot
set tmargin at screen bot+dy
set bmargin at screen bot
rmin=1e1
rmax=1e4
set xlabel 'R / (h^{-1} Mpc)'
set log x
set xrange [rmin:rmax]
set colorbox user origin rig+0.02,bot size 0.03,(top-bot)
plot for [i=1:n] file(i) u 5:6:(2**(i-1)) w l lw 3 lc palette noti

unset multiplot
