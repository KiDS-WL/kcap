reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(print==0){set term aqua}
if(print==1){set term post enh col sol fontfile cmmi; set output 'Cl_redshift.eps'}

cl(z1,z2)=sprintf('data/redshift_%i_%i_cl_full.dat',z1,z2)

ell='{/cmmi10 \140}'

icumulative=1

ellmin=1
ellmax=1e4
set log x
if(print==0){set xlabel 'l'}
if(print==1){set xlabel ell}
set mxtics 10
set xrange [*:*]
set format x '10^{%T}'

set log y
if(print==0){set ylabel 'l(l+1)C_{i,j}(l) / 2{/Symbol p}'}
if(print==1){set ylabel ''.ell.'('.ell.'+1)C_{i,j}('.ell.') / 2{/Symbol p}'}
set format y '10^{%T}'
set mytics 10
#set yrange [1e-9:1e-4]

zmax=1
nz=8
set palette defined (1 'dark-red', 2 'gold')
set cbrange [0:zmax]
set cblabel 'Integration maximum redshift'

set title 'Cross correlation build-up as a function of redshift'

set key top left

if(icumulative==0){
plot for [i=1:nz] cl(i,i+1) u 1:3:(real(i-1)/real(nz-1)) w l lw 3 lc palette noti,\
     'data/cl_full.dat' u 1:3 w l lw 3 lc -1 ti 'Total'
}
if(icumulative==1){
plot for [i=1:nz] cl(0,i) u 1:3:(real(i)/real(nz)) w l lw 3 lc palette noti,\
     'data/cl_full.dat' u 1:3 w l lw 3 lc -1 ti 'Total'
}
