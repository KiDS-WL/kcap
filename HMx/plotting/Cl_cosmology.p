unset multiplot
reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(print==0){set term aqua dashed}
if(print==1){set term post enh col sol fontfile cmmi; set output 'Cl_cosmology.eps'}

file(i)=sprintf('data/cosmology_%d_cl_hm.dat',i)

ell='{/cmmi10 \140}'

pmin=0.7
pmax=0.9
plab='{/Symbol s}_8'
#set palette rgb 34,35,36
set palette defined (1 'blue', 2 'grey', 3 'red')

unset key

set lmargin 10
set rmargin 2

ellmin=1
ellmax=1e4
set log x
set mxtics 10
set xrange [*:*]

set cblabel plab
set cbrange [pmin:pmax]

set multiplot layout 2,1

set xlabel ''
set format x ''

set log y
set mytics 10
if(print==0){set ylabel 'l(l+1)C_{{/Symbol k}{/Symbol k}}(l) / 2{/Symbol p}'}
if(print==1){set ylabel ''.ell.'('.ell.'+1)C_{{/Symbol k}{/Symbol k}}('.ell.') / 2{/Symbol p}'}
set format y '10^{%T}'
#set yrange [1e-8:1e-4]

n=5
plot for [i=1:n] file(i) u 1:3:(pmin+(pmax-pmin)*real(i-1)/real(n-1)) w l lw 3 lc palette dt 1 noti

if(print==0){set xlabel 'l'}
if(print==1){set xlabel ''.ell.''}
set format x '10^{%T}'

unset log y
set yrange [*:*]
set format y
set mytics
if(print==0){set ylabel 'C(l) / C_{mid}(l)'}
if(print==1){set ylabel 'C('.ell.') / C_{mid}('.ell.')'}

plot for [i=1:n] '<paste '.file(i).' '.file((n+1)/2).'' u 1:($3/$6):(pmin+(pmax-pmin)*real(i-1)/real(n-1)) w l lw 3 lc palette dt 1 noti

unset multiplot
