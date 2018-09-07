reset
unset multiplot

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(!exists("print")){print=0}
if(print==0){set term aqua dashed font ',14'; ell='l'}
if(print==1){set term post enh col fontfile cmmi; set output 'Cl.eps'; ell='{/cmmi10 \140}'}

#Files to plot
linear='data/cl_linear.dat'
twohalo='data/cl_2halo.dat'
onehalo='data/cl_1halo.dat'
full='data/cl_full.dat'

#Do linear or not
ilin=0

#Power-law plotting
ifunc=0
ns=0.
#B=1e-10 #Appropriate for kk one-halo term
B=1e-13 #Appropriate for ky one-halo term
#A=1e3
lmax=20
f(l) = l<lmax ? B*(l+1)*l**(1+ns) : 1/0
tits(n)=sprintf('C('.ell.') \~ l^{%1.1f}',n)

#For y axis
#1 - C(l)
#2 - l(l+1)C(l)
#3 - lC(l)/pi - for converting to CMB-lensing potential
if(!exists("icl")){icl=1; print 'Setting icl: ', icl}
if(icl==1){c=2; A=1.; n=0.}
if(icl==2){c=3; A=1.; n=0.}
if(icl==3){c=2; A=1./pi; n=1.}

ellmin=1e0
ellmax=1e4
set log x
set xrange [ellmin:ellmax]
set mxtics 10

set lmargin 10
set rmargin 2

y1=0.98
y2=0.43
y3=0.4
y4=0.1
#bot=0.1
#gap=0.0
#dy=(top-bot-gap)/2.

set multiplot

set tmargin at screen y1
set bmargin at screen y2

set format x ''
set xlabel ''

if(icl==1 || icl==2){set log y; set mytics 10; set format y '10^{%T}'}
if(icl==3){set format y}
if(icl==1){set ylabel 'C_{i,j}('.ell.')'}
if(icl==2){set ylabel ''.ell.'('.ell.'+1)C_{i,j}('.ell.') / 2{/Symbol p}'}
if(icl==3){set ylabel '2'.ell.'C_{i,j}('.ell.') / 2{/Symbol p}'}

#0 - None
#1 - kappa-kappa
#2 - y-y
#3 - CMB-CMB
#4 - kappa-y
#5 - CMB-y
#6 - Comparison to Hill & Spergel (2014) y-phi plot (Fig. 1)
#7 - Comparison to Horowitz & Seljak (2017) y-y plot (Fig. 1)
iaxes=0
if(iaxes==0){set yrange [*:*]}
if(iaxes==1){set yrange [1e-8:1e-4]}
if(iaxes==2){set yrange [1e-15:1e-11]}
if(iaxes==3){set yrange [1e-8:1e-2]}
if(iaxes==4){set yrange [1e-12:1e-8]}
if(iaxes==5){set yrange [1e-12:1e-7]}
if(iaxes==6){set yrange [0.:4e-11]}
if(iaxes==7){set yrange [1e-14:1e-11]}

set key top left

plot twohalo u 1:(column(c)*A*$1**n) w l lw 3 lc 'red' dt 2 ti '2-halo term',\
     onehalo u 1:(column(c)*A*$1**n) w l lw 3 lc 'red' dt 3 ti '1-halo term',\
     full    u 1:(column(c)*A*$1**n) w l lw 3 lc 'red' dt 1 ti 'Full'#,\
     #f(x)                            w l lw 2 lc -1    dt 2 ti tits(ns),\
     #linear  u 1:(column(c)*A*$1**n) w l lw 3 lc -1    dt 1 ti 'Linear'

if(ifunc==1){replot f(x) w l lw 2 lc -1 dt 2 ti tits(ns)}
if(ilin==1) {replot linear  u 1:(column(c)*A*$1**n) w l lw 3 lc -1 dt 1 ti 'Linear'}

set tmargin at screen y3
set bmargin at screen y4

set format x '10^{%T}'
if(print==0) set xlabel 'l'
if(print==1) set xlabel ''.ell.''

unset log y
set yrange [0:1]
set format y
set mytics
if(print==0) set ylabel 'C_{n,halo}(l) / C_{full}(l)'
if(print==1) set ylabel 'C_{n,halo}('.ell.') / C_{full}('.ell.')'

plot '<paste '.twohalo.' '.full.'' u 1:($3/$6) w l lw 3 lc 'red' dt 2 noti,\
     '<paste '.onehalo.' '.full.'' u 1:($3/$6) w l lw 3 lc 'red' dt 3 noti

unset multiplot
