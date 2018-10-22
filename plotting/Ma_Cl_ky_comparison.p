reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(print==0) set term aqua dashed
if(print==1) set term post enh col fontfile cmmi; set output 'Ma_comparison.eps'

#Non-standard font
if(print==0) ell='l'
if(print==1) ell='{/cmmi10 \140}'

#if(print==0) set xlabel 'l'
#if(print==1) set xlabel ''.ell.''
ellmin=1
ellmax=1e4
set xlabel ''.ell.''
set log x
set format x '10^{%T}'
set xrange [ellmin:ellmax]
set mxtics 10

#if(print==0) set ylabel 'l(l+1)C_{y{/Symbol k}}(l) / 2{/Symbol p}'
#if(print==1) set ylabel ''.ell.'('.ell.'+1)C_{y{/Symbol k}}('.ell.') / 2{/Symbol p}'
set ylabel ''.ell.'('.ell.'+1) C_{y{/Symbol k}}('.ell.') / 2{/Symbol p}'
set log y
set yrange [1e-13:1e-8]
set format y '10^{%T}'
set mytics 10

HM_Planck='RCSLenS/Cell_ky_planck_RCS.dat'
HM_WMAP='RCSLenS/Cell_ky_wmap_RCS.dat'
HM_Mead='data/cl_hm.dat'

set key top left

set title 'Comparison with Ma et al. (2015) y-{/Symbol k} cross correlation'

#power-law stuff
fac=1e-13
ns=0.
lmax=20
f(l) = l<lmax ? fac*(l+1)*l**(1+ns) : 1/0

tits(n)=sprintf('C('.ell.') ~ l^{%1.1f}',n)

plot HM_Planck u 1:($2) w l dt 1 lc 1 lw 3 ti 'Halo model - Ma (2015) - Planck',\
     HM_Planck u 1:($3) w l dt 3 lc 1 lw 3 noti,\
     HM_Planck u 1:($4) w l dt 4 lc 1 lw 3 noti,\
     HM_WMAP   u 1:($2) w l dt 1 lc 2 lw 3 ti 'Halo model - Ma (2015) - WMAP',\
     HM_WMAP   u 1:($3) w l dt 3 lc 2 lw 3 noti,\
     HM_WMAP   u 1:($4) w l dt 4 lc 2 lw 3 noti,\
     'data/cl_hm.dat' u 1:($3) w l lw 3 dt 1 lc 3 ti 'Halo model - Mead',\
     'data/cl_1h.dat' u 1:($3) w l lw 3 dt 3 lc 3 noti,\
     'data/cl_2h.dat' u 1:($3) w l lw 3 dt 4 lc 3 noti,\
      f(x)                        w l lw 3 dt 2 lc -1 ti tits(ns)
