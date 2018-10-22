reset
unset multiplot

set term aqua size 1000,1000

# Initial white space
print ''

# Data file
if(!exists('chain')){chain='fitting/test_chain.dat'}
print 'Plotting data file: chain: ', chain

chains="'fitting/AGN_TUNED_nu0_z0.0_n50000_m13_chain.dat' 'fitting/AGN_7p6_nu0_z0.0_n50000_m13_chain.dat' 'fitting/AGN_8p0_nu0_z0.0_n50000_m13_chain.dat'"

if(!exists('triangle')){triangle=1}
print 'triangle: ', triangle

# Chain pruning
if(!exists('m')){m=0}
print 'Skipping entries: m: ', m
print ''

# Binning function
min=0.
n=1000
width = 1./real(n) # binwidth; evaluates to 1.0
bin(x) = min+width*(floor((x-min)/width)+0.5)

if(!exists('im')){im=15}
print 'im: 11 - Everything'
print 'im: 12 - Gas'
print 'im: 13 - Stars'
print 'im: 14 - Gas and stars'
print 'im: 15 - Matter'
print 'im: ', im
print ''

# Everything
if(im==11) {
np=8
labels="'log_{10}{/Symbol a}' 'log_{10}{/Symbol e}' '{/Symbol G}-1' 'log_{10} M_0' 'log_{10} A_*' 'log_{10} T_{w}' 'log_{10} c_*' 'log_{10} f_c'"
}

# Gas
if(im==12) {
np=6
labels="'log_{10}{/Symbol e}' '{/Symbol G}-1' 'log_{10} M_0' 'log_{10} A_*' 'log_{10} f_c' '{/Symbol G}^p'"
}

# Stars
if(im==13) {
np=5
labels="'log_{10}{A_*}' 'log_{10}{c_*}' 'c_*^p' 'log_{10}M_*' '{/Symbol s}_*'"
}

# Gas and stars
if(im==14) {

}

# Matter
if(im==15) {
np=3
labels="'log_{10}{/Symbol e}' '{/Symbol G}-1' 'log_{10} M_0' 'log_{10} A_*' 'log_{10} c_*' 'log_{10} f_c' '{/Symbol G}^p' 'c_*^p' 'log_{10}M_*' '{/Symbol s}_*'"
}

# Number of parameters
print 'Number of parameters: np: ', np
print ''

# Write parameters to screen
print 'Parameters to plot:'
do for [i=1:np] {print 'Parameter: ', i, ' ', word(labels,i)}
print ''

# Plot boundary stuff
top=0.98
bot=0.10
lef=0.10
rig=0.98
nx=np
ny=np
dx=(rig-lef)/ny
dy=(top-bot)/nx

unset colorbox
set palette grey negative

# Now actually do the plotting
set multiplot

#do for [i=1:np-1]{
do for [i=1:np]{

do for [j=i:np]{

set lmargin at screen lef+(i-1)*dx
set rmargin at screen lef+(i-0)*dx

set tmargin at screen top-(j-1)*dy
set bmargin at screen top-(j-0)*dy

if(i==1)  {set format y; set ylabel word(labels,j); set ytics} else {set format y ''; set ylabel ''}
if(j==np) {set format x; set xlabel word(labels,i); set ytics} else {set format x ''; set xlabel ''}
if(i==j)  {set format y ''; set ylabel ''; unset ytics}
if(i==j && triangle==0) {set xlabel word(labels,i); set format x}

if(i==j){
plot for [k=1:words(chains)] word(chains,k) every ::m u (bin(column(i+1))):(1.0) lc k smooth freq with boxes noti
}

#if(i!=j){
#plot chain every ::m u (column(i+1)):(column(j+1)) w d lc -1 noti
#}

if(triangle==1 && i!=j){
#set pm3d at b    # draw on bottom, not as 3d surface
#set view map     # don't do a 3-d looking plot
#set dgrid 30,30  # grid of 100x100 pixels
#splot chain u (column(i+1)):(column(j+1)):1 w pm3d noti
plot for [k=1:words(chains)] word(chains,k) every ::m u (column(i+1)):(column(j+1)) w d lc k noti
}

}
}

unset multiplot
