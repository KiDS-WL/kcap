reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(!exists('print')){print=0}
if(print==0){set term aqua dashed; ell='l'}
if(print==1){set term post enh col dashed fontfile cmmi; set output 'triad.eps'; ell='{/cmmi10 \140}'}

#File name and location functions
sim_file(x,y,mod)=sprintf('/Users/Mead/Physics/people/Tilman/Cl/Cl_%s-%s_%s.txt', x, y, mod)
hm_file(x,y)=sprintf('data/triad_Cl_%s-%s.dat',x,y)

#Neutrino simulations (BAHAMAS)
#sims="'0.00_eV' '0.06_eV' '0.12_eV' '0.24_eV' '0.48_eV'"
#names="'0.00 eV' '0.06 eV' '0.12 eV' '0.24 eV' '0.48 eV'"

#Feedback simulations (cosmo-OWLS)
#sims="'REF' 'NOCOOL' 'AGN_8.0' 'AGN_8.5' 'AGN_8.7'"
#names="'REF' 'NOCOOL' 'AGN 8.0' 'AGN 8.5' 'AGN 8.7'"
sims="'AGN_8.0' 'AGN_8.5' 'AGN_8.7'"
names="'AGN 8.0' 'AGN 8.5' 'AGN 8.7'"

#Types of fields
things="'CMB' 'y' 'gal'"

#Symbols for fields for legend (e.g., phi, y, gamma)
symbols="'{/Symbol f}' 'y' '{/Symbol g}'"

#x axis
ellmin=90.
ellmax=3500.
set log x
set xlabel ''.ell.''
set xrange [ellmin:ellmax]

#y axis
clmin=1e-14
clmax=2e-12
set log y
#set ylabel 'C_{i,j}('.ell.')'
#set yrange [clmin:clmax]
set yrange [*:*]

#Distance and function to shift simulation points for clarity
#dmax=1.1
disp(i,n)=1.+0.03*real(i-1)/real(n-1)

#Factor to shift down CMB-gal curves
f=5e-4

#1 - C(l)
#2 - l(l+1)C(l)/2pi
icl=2

#Number of simulation curves to make
n=words(sims)

if(icl==1){p1=0; p2=0; p3=0; ylab='C_{i,j}('.ell.')'                                 ; c=2; clmin=1e-14; clmax=2e-12; set key top right}
if(icl==2){p1=1; p2=1; p3=1; ylab=''.ell.'('.ell.'+1)C_{i,j}('.ell.') / 2{/Symbol p}'; c=3; clmin=1e-10; clmax=1e-7; set key bottom right}

set ylabel ylab
set yrange [clmin:clmax]
#Do the plot
plot hm_file(word(things,3),word(things,1)) u 1:(f*column(c)) w l lw 3 lc -1 dt 1 ti ''.word(symbols,3).'-'.word(symbols,1).'',\
     hm_file(word(things,1),word(things,2)) u 1:(column(c))   w l lw 3 lc -1 dt 2 ti ''.word(symbols,1).'-'.word(symbols,2).'',\
     hm_file(word(things,2),word(things,3)) u 1:(column(c))   w l lw 3 lc -1 dt 3 ti ''.word(symbols,2).'-'.word(symbols,3).'',\
     for [i=1:n] sim_file(word(things,3),word(things,1),word(sims,i)) u ($1*disp(i,n)):(f*($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(f*($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3) w errorbars lc i pt 7 noti,\
     for [i=1:n] sim_file(word(things,1),word(things,2),word(sims,i)) u ($1*disp(i,n)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3)     w errorbars lc i pt 7 ti word(names,i),\
     for [i=1:n] sim_file(word(things,2),word(things,3),word(sims,i)) u ($1*disp(i,n)):(($1**p1)*(($1+1)**p2)*$2/(2.*pi)**p3):(($1**p1)*(($1+1)**p2)*$3/(2.*pi)**p3)     w errorbars lc i pt 7 noti

#l(l+1)C(l)
#if(icl==2){
#Set the key
#set key bottom right
#y axis
#clmin=1e-10
#clmax=1e-7
#set ylabel ''.ell.'('.ell.'+1)C('.ell.') / 2{/Symbol p}'
#set yrange [clmin:clmax]
#Do the plot
#plot hm_file(word(things,3),word(things,1)) u 1:(f*$3) w l lw 3 lc -1 dt 1 ti 'gal-CMB',\
#     hm_file(word(things,1),word(things,2)) u 1:3      w l lw 3 lc -1 dt 2 ti 'CMB-y',\
#     hm_file(word(things,2),word(things,3)) u 1:3      w l lw 3 lc -1 dt 3 ti 'y-gal',\#
#     for [i=1:n] sim_file(word(things,3),word(things,1),word(sims,i)) u w errorbars lc i pt i noti,\
#     for [i=1:n] sim_file(word(things,1),word(things,2),word(sims,i)) u w errorbars lc i pt i ti word(names,i),\
#     for [i=1:n] sim_file(word(things,2),word(things,3),word(sims,i)) u    w errorbars lc i pt i noti
#}
