reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(!exists('print')){print=0}
if(print==0){set term aqua dashed; ell='l'}
if(print==1){set term post enh col fontfile cmmi font ',12'; set output 'cosmoSIS_triad_comparisons.eps'; ell='{/cmmi10 \140}'}

hmx(f1,f2)=sprintf('data/triad_Cl_%s-%s.dat',f1,f2)

cs(f1,f2,thing)=sprintf('/Users/Mead/Physics/people/Tilman/cosmoSIS/%s_%s_cl/%s.txt',f1,f2,thing)

lmin=50
lmax=5000
set log x
xlab=''.ell.''
set xrange [lmin:lmax]

lef=0.1
rig=0.98
xmid=(lef+rig)/2.
dx=0.03

top=0.98
bot=0.10
ymid=(top+bot)/2.

ylab=''.ell.'('.ell.'+1)C('.ell.') / 2{/Symbol p}'

set multiplot layout 2,2

## Shear - CMB lensing (top left) ##

set tmargin at screen top
set bmargin at screen ymid
set lmargin at screen lef
set rmargin at screen xmid-dx

set format x ''
set xlabel ''

clmin=5e-6
clmax=2e-4
set log y
set ylabel ylab
set format y '10^{%T}'
set yrange [clmin:clmax]

set key top left

plot hmx('gal_z0.1-0.9','CMB') u 1:3 w l lw 2 lc 1 ti 'HMx: {/Symbol g} (z = 0.1->0.9)-{/Symbol k}',\
     hmx('gal_z0.1-0.5','CMB') u 1:3 w l lw 2 lc 2 ti 'HMx: {/Symbol g} (z = 0.1->0.5)-{/Symbol k}',\
     hmx('gal_z0.5-0.9','CMB') u 1:3 w l lw 2 lc 3 ti 'HMx: {/Symbol g} (z = 0.5->0.9)-{/Symbol k}',\
     '<paste '.cs('shear','cmbkappa','ell').' '.cs('shear','cmbkappa','bin_1_1').'' u 1:($1*($1+1)*$2/(2*pi)) w l lw 2 lc 1 dt 2 ti 'cosmoSIS: {/Symbol g} (z = 0.1->0.9)-{/Symbol k}',\
     '<paste '.cs('shear','cmbkappa','ell').' '.cs('shear','cmbkappa','bin_2_1').'' u 1:($1*($1+1)*$2/(2*pi)) w l lw 2 lc 2 dt 2 ti 'cosmoSIS: {/Symbol g} (z = 0.1->0.5)-{/Symbol k}',\
     '<paste '.cs('shear','cmbkappa','ell').' '.cs('shear','cmbkappa','bin_3_1').'' u 1:($1*($1+1)*$2/(2*pi)) w l lw 2 lc 3 dt 2 ti 'cosmoSIS: {/Symbol g} (z = 0.5->0.9)-{/Symbol k}'

## ##

## Shear - y (top right) ##

set tmargin at screen top
set bmargin at screen ymid
set lmargin at screen xmid+dx
set rmargin at screen rig

set xlabel xlab
set format x

clmin=1e-10
clmax=1e-8
set log y
set ylabel ''
set format y '10^{%T}'
set yrange [clmin:clmax]

set key top left

plot hmx('y','gal_z0.1-0.9') u 1:3 w l lw 2 lc 1 ti 'HMx: {/Symbol g} (z = 0.1->0.9)-y',\
     hmx('y','gal_z0.1-0.5') u 1:3 w l lw 2 lc 2 ti 'HMx: {/Symbol g} (z = 0.1->0.5)-y',\
     hmx('y','gal_z0.5-0.9') u 1:3 w l lw 2 lc 3 ti 'HMx: {/Symbol g} (z = 0.5->0.9)-y',\
     '<paste '.cs('shear','y','ell').' '.cs('shear','y','bin_1_1').'' u 1:($1*($1+1)*$2/(2*pi)) w l lw 2 lc 1 dt 2 ti 'cosmoSIS: {/Symbol g} (z = 0.1->0.9)-y',\
     '<paste '.cs('shear','y','ell').' '.cs('shear','y','bin_2_1').'' u 1:($1*($1+1)*$2/(2*pi)) w l lw 2 lc 2 dt 2 ti 'cosmoSIS: {/Symbol g} (z = 0.1->0.5)-y',\
     '<paste '.cs('shear','y','ell').' '.cs('shear','y','bin_3_1').'' u 1:($1*($1+1)*$2/(2*pi)) w l lw 2 lc 3 dt 2 ti 'cosmoSIS: {/Symbol g} (z = 0.5->0.9)-y'

## ##

## CMB lensing-y (bottom left) ##

set tmargin at screen ymid
set bmargin at screen bot
set lmargin at screen lef
set rmargin at screen xmid-dx

set xlabel xlab
set format x

clmin=3e-10
clmax=3e-8
set log y
set ylabel ylab
set format y '10^{%T}'
set yrange [clmin:clmax]

set key top left

f=1.

plot hmx('CMB','y') u 1:3 w l lw 2 lc 1 ti 'HMx: {/Symbol k}-y',\
     '<paste '.cs('cmbkappa','y','ell').' '.cs('cmbkappa','y','bin_1_1').'' u 1:(f*$1*($1+1)*$2/(2*pi)) w l lw 2 lc 1 dt 2 ti 'cosmoSIS: {/Symbol k}-y'

##

unset multiplot
