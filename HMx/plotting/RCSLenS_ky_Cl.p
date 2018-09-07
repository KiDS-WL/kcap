reset

cmmi='/Users/Mead/Fonts/cmmi10.pfb'

if(print==0) set term aqua dashed
if(print==1) set term post enh col fontfile cmmi; set output 'RCSLenS_ky_Cl.eps'

#Non-standard font
ell='{/cmmi10 \140}'

#Simulation C(l) are multiplied by this number for reasons unknown
simfac=1e13

#Sigmas for the Planck beam
sigma1 = 0.00742*10./60.
sigma2 = 0.00742*9.5/60.

kernel1(l) = exp(-0.5*l*(l+1.)*sigma1**2.)
kernel2(l) = exp(-0.5*l*(l+1.)*sigma2**2.)# same as: Planck_smoothing = np.exp(-sigma_Planck**2*ell_ky**2/2.0
kernel(l) = kernel1(l)*kernel2(l)

if(print==0) set xlabel 'l'
if(print==1) set xlabel ''.ell.''
set xrange [0:2000]

if(print==0) set ylabel 'l(l+1)C_{y{/Symbol k}}(l) / 2{/Symbol p}'
if(print==1) set ylabel ''.ell.'('.ell.'+1)C_{y{/Symbol k}}('.ell.') / 2{/Symbol p}'
set yrange [0:2e-9]

RCSLenS='RCSLenS/Cell_y-kappa_RCS-2.txt'
sims='RCSLenS/Cell_y-kappa_sims.txt'
HM_Planck='RCSLenS/Cell_ky_planck_RCS.dat'
HM_WMAP='RCSLenS/Cell_ky_wmap_RCS.dat'
HM_Mead='data/cl_full.dat'

set key top right

set title 'RCSLenS y-{/Symbol k} cross correlation'

plot 0 w l lt -1 noti,\
     sims u 1:($1*($1+1)*$2/((2.*pi)*simfac)) w l lw 3 dt 2 lc rgb '#BB000000' ti 'AGN - WMAP',\
     sims u 1:($1*($1+1)*$3/((2.*pi)*simfac)) w l lw 3 dt 2 lc 4 ti 'AGN - Planck',\
     sims u 1:($1*($1+1)*$4/((2.*pi)*simfac)) w l lw 3 dt 2 lc 5 ti 'NOCOOL - WMAP',\
     sims u 1:($1*($1+1)*$5/((2.*pi)*simfac)) w l lw 3 dt 2 lc 6 ti 'REF - WMAP',\
     sims u 1:($1*($1+1)*$6/((2.*pi)*simfac)) w l lw 3 dt 2 lc 7 ti 'AGN 8.5 - WMAP',\
     sims u 1:($1*($1+1)*$7/((2.*pi)*simfac)) w l lw 3 dt 2 lc 8 ti 'AGN 8.7 - WMAP',\
     HM_Planck u 1:($2*kernel($1)) w l dt 1 lc 1 lw 3 ti 'Halo model - Ma (2015) - Planck',\
     HM_WMAP   u 1:($2*kernel($1)) w l dt 1 lc 2 lw 3 ti 'Halo model - Ma (2015) - WMAP',\
     HM_WMAP   u 1:($3*kernel($1)) w l dt 3 lc 2 lw 3 noti,\
     HM_WMAP   u 1:($4*kernel($1)) w l dt 4 lc 2 lw 3 noti,\
     'data/cl_full.dat'  u 1:($3*kernel($1)) w l lw 3 dt 1 lc 3 ti 'Halo model - Mead',\
     'data/cl_1halo.dat' u 1:($3*kernel($1)) w l lw 3 dt 3 lc 3 noti,\
     'data/cl_2halo.dat' u 1:($3*kernel($1)) w l lw 3 dt 4 lc 3 noti,\
     RCSLenS u 1:($1*($1+1)*$2/(2.*pi)):($1*($1+1)*$3/(2.*pi)) w errorbars ps 1 pt 7 lc rgb 'black' ti 'RCSLenS data'
