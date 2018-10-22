reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'
cmsy='/Users/Mead/Fonts/cmsy10.pfb'

if(print==0){set term aqua; ell='l'; msun='M'}
if(print==1){set term post enh col sol fontfile cmmi fontfile cmsy; set output 'Cl_mass.eps'; ell='{/cmmi10 \140}'; msun='{/cmsy10 \014}'}

cl(m1,m2)=sprintf('data/mass_%i_%i_cl_hm.dat',m1,m2)

icumulative=1

ellmin=1
ellmax=1e4
set log x
set xlabel ''.ell.''
set mxtics 10
set xrange [*:*]
set format x '10^{%T}'

set log y
set ylabel ''.ell.'('.ell.'+1)C_{i,j}('.ell.') / 2{/Symbol p}'
set format y '10^{%T}'
set mytics 10

mmin=1e11
mmax=1e16
nm=6
set palette defined (1 'pink', 2 'dark-red')
set log cb
set cbrange [mmin:mmax]
set format cb '10^{%T}'
set cblabel 'Integration maximum mass [M_'.msun.' / h]'

set key top left

set title 'Cross correlation build-up as a function of halo mass'

if(icumulative==0){
plot for [i=10:15] cl(i,i+1) u 1:3 w l lw 3 lc i ti tits(i,i+1),\
     'data/cl_hm.dat' u 1:3 w l lw 3 lc -1 ti 'Total'
}

if(icumulative==1){
plot for [i=1:nm] cl(7,i+10) u 1:3:(exp(log(mmin)+(log(mmax/mmin))*real(i-1)/real(nm-1))) w l lw 3 lc palette noti#,\
     'data/cl_hm.dat' u 1:3 w l lw 3 lc -1 ti 'Total'
}

