unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed}
if(print==1) {set term post enh col dashed dl .5 font ',10'}

print ''

# Simulations to compare against
# 1 - cosmo-OWLS
# 2 - BAHAMAS
# 3 - Generic hydro comparison
if(!exists('icomp')){icomp=2}
print 'icomp = 1: Compare to cosmo-OWLS (NO LONGER SUPPORTED)'
print 'icomp = 2: Compare to BAHAMAS'
print 'icomp = 3: Generic hydro, no comparison simulation'
print 'icomp = '.icomp.''
print ''

if(icomp==1){print 'Twat; icomp=1 does not work'; exit}

# if(icomp==1){sims='cosmo-OWLS'; Om_m=0.272; Om_b=0.0455}
if(icomp==2){sims='BAHAMAS'}
if(icomp==3){sims=''}

# cosmological parameters (only used for plotting Om_b/Om_m lines)
Om_m=0.2793
Om_b=0.0463
Om_c=Om_m-Om_b
print 'Omega_m: ', Om_m
print 'Omega_b: ', Om_b
print 'Omega_c: ', Om_c
print ''

# Redshift
if(!exists('z')){z=0.0}
if(z==0.0){snap='snap32'}
if(z==0.5){snap='snap28'}
if(z==1.0){snap='snap26'}
if(z==2.0){snap='snap22'}
print 'Redshift: z = ', z
print ''

# Plot to make
if(!exists('iplot')){iplot=1}
print 'iplot = 1: Power spectrum plot'
print 'iplot = 2: Power spectrum ratio plot'
print 'iplot = 3: Power spectrum suppression plot'
print 'iplot = 4: Power spectrum residual plot'
print 'iplot = 5: Power spectrum components'
print 'iplot = 6: Power spectrum of electron pressure'
print 'iplot = 7: Power spectrum with k^1.5 units'
print 'iplot = '.iplot.''
print ''

# File names - cosmo-OWLS
#if(icomp==1){
#data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',sim,type1,type2)
#hmpk(sim,i,j)=sprintf('cosmo-OWLS/power_%s_%i%i.dat',sim,i,j)
#print 'Error: cosmo-OWLS comparison needs to be updated'
#exit
#}

if(!exists('mesh')){mesh=1024}
print 'Mesh size: *mesh*: ', mesh
print ''

# Simulation data files
plot_title_name_z(sim,z)=sprintf('BAHAMAS comarison of %s at z = %1.1f', sim, z)
plot_title_z(z)=sprintf('BAHAMAS comarison at z = %1.1f', z)
data(sim,mesh,s,type1,type2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M%d/%s_nu0_L400N1024_WMAP9_%s_%s_%s_power.dat',mesh,sim,s,type1,type2)
data_dmonly='DMONLY_2fluid'

# File names - BAHAMAS
if(icomp==2){
hmpk(sim,z,i,j)=sprintf('BAHAMAS/power_%s_z%1.1f_%i%i.dat',sim,z,i,j)
hmpk_dmonly=hmpk('DMONLY',z,0,0)
}

# File names - generic hydro
if(icomp==3){
hmpk(sim,z,i,j)=sprintf('hydro/power_z%1.1f_%i%i.dat',z,i,j)
hmdm(z)=sprintf('hydro/power_z%1.1f.dat',z)
hmpk_dmonly=hmdm(z)
}

# Columns for simulation power
c=2
s=3
L=5

# Columns for halo-model power
d=5
M=5

#cosmo-OWLS simulation names
#if(icomp==1){
#hm_names="'DMONLY' 'REF' 'NOCOOL' 'AGN' 'AGN8p5' 'AGN8p7'"
#owls_names="'DMONLY' 'REF' 'NOCOOL_UVB' 'AGN' 'AGN_Theat_8p5' 'AGN_Theat_8p7'"
#}

# BAHAMAS simulation names
#if(icomp==2){hmpk_names="'DMONLY' 'AGN-lo' 'AGN-hi' 'AGN'"}
if(icomp==2){hmpk_names="'DMONLY' 'AGN-lo' 'AGN' 'AGN-hi'"}
if(icomp==3){hmpk_names="''"}
#data_names="'DMONLY_2fluid' 'AGN_7p6' 'AGN_8p0' 'AGN_TUNED'"
data_names="'DMONLY_2fluid' 'AGN_7p6' 'AGN_TUNED' 'AGN_8p0'"

# Set the comparison model
if(!exists('nsim')){nsim=3}
hmpk_name=word(hmpk_names,nsim)
print 'Variable: *nsim* '.nsim.''
print 'Halo model power file name: '.hmpk_name.''
data_name=word(data_names,nsim)
print 'Simuation power file: '.data_name.''
print ''

# All different fields for power spectra
fld0='all'
fld1='dm'
fld2='gas'
fld3='stars'
fld6='epressure'

print 'Example simulation file: ', data(data_name,mesh,snap,fld0,fld0)
print 'Example halo-model file: ', hmpk(hmpk_name,z,0,0)
print ''

# Factor to adject the (dimensionful) pressure spectrum
f=1e2
print 'Pressure spectra multiplied by: ', f
print ''

# Set colours
col1=1 # All matter
col2=2 # Dark matter
col3=3 # Gas
col4=4 # Stars
col5=7 # Matter-pressure
col6=6 # Pressure

# k range
if(icomp==1 || icomp==2){kmin=1e-2; kmax=1e1}
if(icomp==3){kmin=1e-2; kmax=1e1}
klab='k / h Mpc^{-1}'
set xlabel klab
set format x
set log x
set xrange [kmin:kmax]

# Delta^2(k) range
pmin=1e-7
pmax=1e3
set log y
set yrange [pmin:pmax]
set format y '10^{%T}'
set ylabel '{/Symbol D}_{i,j}^2(k)'
set mytics 10

# Set the overall plot titles
set title plot_title_name_z(data_name,z)

if(iplot==1){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
set output outfile(data_name,z)
}

#set title 'Comparison of '.sims.' '.hm_name.' simulation to halo-model predictions'

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     data(data_name,mesh,snap,fld0,fld0) u 1:(column(c)-column(s)):5 w e pt 7 lc col1 noti,\
     data(data_name,mesh,snap,fld1,fld1) u 1:(column(c)-column(s)):5 w e pt 7 lc col2 noti,\
     data(data_name,mesh,snap,fld2,fld2) u 1:(column(c)-column(s)):5 w e pt 7 lc col3 noti,\
     data(data_name,mesh,snap,fld3,fld3) u 1:(column(c)-column(s)):5 w e pt 7 lc col4 noti,\
     data(data_name,mesh,snap,fld0,fld1) u 1:(column(c)-column(s)):5 w e pt 6 lc col2 noti,\
     data(data_name,mesh,snap,fld0,fld2) u 1:(column(c)-column(s)):5 w e pt 6 lc col3 noti,\
     data(data_name,mesh,snap,fld0,fld3) u 1:(column(c)-column(s)):5 w e pt 6 lc col4 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     hmpk(hmpk_name,z,1,1) u 1:(column(d)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(hmpk_name,z,2,2) u 1:(column(d)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     hmpk(hmpk_name,z,3,3) u 1:(column(d)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     hmpk(hmpk_name,z,0,1) u 1:(column(d)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(hmpk_name,z,0,2) u 1:(column(d)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(hmpk_name,z,0,3) u 1:(column(d)) w l lw 3 dt 2 lc col4 noti

}

if(iplot==7){

pow=1.5

# Delta^2(k)/k^1.5 range
pmin=1e-4
pmax=1e2
set log y
set yrange [pmin:pmax]
set format y '10^{%T}'
set ylabel '{/Symbol D}_{i,j}^2(k) / [k / h^{-1} Mpc]^{1.5}'
set mytics 10

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
set output outfile(data_name,z)
}

#set title 'Comparison of '.sims.' '.hm_name.' simulation to halo-model predictions'

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     data(data_name,mesh,snap,fld0,fld0) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col1 noti,\
     data(data_name,mesh,snap,fld1,fld1) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col2 noti,\
     data(data_name,mesh,snap,fld2,fld2) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col3 noti,\
     data(data_name,mesh,snap,fld3,fld3) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col4 noti,\
     data(data_name,mesh,snap,fld0,fld1) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 6 lc col2 noti,\
     data(data_name,mesh,snap,fld0,fld2) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 6 lc col3 noti,\
     data(data_name,mesh,snap,fld0,fld3) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 6 lc col4 noti,\
     hmpk(hmpk_name,z,0,0) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     hmpk(hmpk_name,z,1,1) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(hmpk_name,z,2,2) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     hmpk(hmpk_name,z,3,3) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     hmpk(hmpk_name,z,0,1) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(hmpk_name,z,0,2) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(hmpk_name,z,0,3) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 2 lc col4 noti

}

if(iplot==6){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_pressure.eps',name,z)
set output outfile(data_name,z)
}

set multiplot layout 2,1

unset xlabel
set format x ''

set key top left

plot data(data_name,mesh,snap,fld0,fld0) u 1:(column(c)-column(s)):5              w e pt 7 lc col1 noti,\
     data(data_name,mesh,snap,fld0,fld6) u 1:(f*(column(c)-column(s))):(f*$5)     w e pt 6 lc col5 noti,\
     data(data_name,mesh,snap,fld6,fld6) u 1:(f*f*(column(c)-column(s))):(f*f*$5) w e pt 7 lc col6 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d))       w l lw 3 dt 1 lc col1 ti 'Matter-Matter',\
     hmpk(hmpk_name,z,0,0) u 1:(column(d-2))     w l lw 3 dt 2 lc col1 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d-1))     w l lw 3 dt 3 lc col1 noti,\
     hmpk(hmpk_name,z,0,6) u 1:(f*column(d))     w l lw 3 dt 1 lc col5 ti 'Matter-Pressure',\
     hmpk(hmpk_name,z,0,6) u 1:(f*column(d-2))   w l lw 3 dt 2 lc col5 noti,\
     hmpk(hmpk_name,z,0,6) u 1:(f*column(d-1))   w l lw 3 dt 3 lc col5 noti,\
     hmpk(hmpk_name,z,6,6) u 1:(f*f*column(d))   w l lw 3 dt 1 lc col6 ti 'Pressure-Pressure',\
     hmpk(hmpk_name,z,6,6) u 1:(f*f*column(d-2)) w l lw 3 dt 2 lc col6 noti,\
     hmpk(hmpk_name,z,6,6) u 1:(f*f*column(d-1)) w l lw 3 dt 3 lc col6 noti

unset title

set xlabel klab
set format x

# Set the y axis for P(k) / P_{DMONLY}(k)
rmin=1e-5
rmax=2
set yrange [rmin:rmax]
set ylabel 'P(k) / P_{DMONLY}(k)'

plot 1 w l lt -1 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld0).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L)))     w p pt 7 lc col1 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld6).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:(f*(column(c)-column(s))/(column(c+L)-column(s+L)))   w p pt 6 lc col5 noti,\
     '<paste '.data(data_name,mesh,snap,fld6,fld6).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:(f*f*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col6 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M))     w l lw 3 dt 1 lc col1 noti '{/Symbol d}{/Symbol d}',\
     '<paste '.hmpk(hmpk_name,z,0,6).' '.hmpk_dmonly.'' u 1:(f*column(d)/column(d+M))   w l lw 3 dt 1 lc col5 noti '{/Symbol d}p',\
     '<paste '.hmpk(hmpk_name,z,6,6).' '.hmpk_dmonly.'' u 1:(f*f*column(d)/column(d+M)) w l lw 3 dt 1 lc col6 noti 'pp'

unset multiplot

}

if(iplot==2){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_ratio.eps',name,z)
set output outfile(data_name,z)
}

set multiplot layout 1,2

# Set the ratio axis
rmin=2e-3
rmax=1.5
set log y
set yrange [rmin:rmax]
set ylabel 'P_{OWL}/P_{DMONLY}'
set mytics 10

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     Om_b/Om_m      w l lc -1 dt 2 noti,\
     Om_c/Om_m      w l lc -1 dt 2 noti,\
     (Om_b/Om_m)**2 w l lc -1 dt 2 noti,\
     (Om_c/Om_m)**2 w l lc -1 dt 2 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld0).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col1 noti,\
     '<paste '.data(data_name,mesh,snap,fld1,fld1).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col2 noti,\
     '<paste '.data(data_name,mesh,snap,fld2,fld2).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col3 noti,\
     '<paste '.data(data_name,mesh,snap,fld3,fld3).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col4 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld0).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col1 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld1).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col2 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld2).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col3 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld3).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col4 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'All matter',\
     '<paste '.hmpk(hmpk_name,z,1,1).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     '<paste '.hmpk(hmpk_name,z,2,2).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col3 ti 'Gas',\
     '<paste '.hmpk(hmpk_name,z,3,3).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col4 ti 'Stars',\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col1 noti,\
     '<paste '.hmpk(hmpk_name,z,0,1).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col2 noti,\
     '<paste '.hmpk(hmpk_name,z,0,2).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col3 noti,\
     '<paste '.hmpk(hmpk_name,z,0,3).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col4 noti

rmin=1e-6
rmax=2.
set yrange [rmin:rmax]
set format y '10^{%T}'

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 noti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 noti 'Cross with matter',\
     '<paste '.data(data_name,mesh,snap,fld0,fld0).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col1 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld6).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((f*column(c))/(column(c+L)-column(s+L)))         w p pt 6 lc col5 noti,\
     '<paste '.data(data_name,mesh,snap,fld6,fld6).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((f*f*column(c))/(column(c+L)-column(s+L)))       w p pt 7 lc col6 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M))     w l lw 3 dt 1 lc col1 ti 'Matter-Matter',\
     '<paste '.hmpk(hmpk_name,z,0,6).' '.hmpk_dmonly.'' u 1:(f*column(d)/column(d+M))   w l lw 3 dt 2 lc col5 ti 'Matter-Pressure',\
     '<paste '.hmpk(hmpk_name,z,6,6).' '.hmpk_dmonly.'' u 1:(f*f*column(d)/column(d+M)) w l lw 3 dt 1 lc col6 ti 'Pressure-Pressure'

unset multiplot

}

if(iplot==3){

set title plot_title_z(z)

if(print==1){
outfile(z)=sprintf('suppression_z%1.1f.eps',z)
set output outfile(z)
}

rmin=0.75
rmax=1.1
unset log y
set format y
set yrange [rmin:rmax]
set ylabel 'P(k) / P_{DMONLY}(k)'

cols="'black' 'orange' 'orange-red' 'red'"

if(icomp==1 || icomp==2){
plot 1 w l lt -1 noti,\
     for [i=1:words(data_names)] '<paste '.data(word(data_names,i),mesh,snap,fld0,fld0).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 dt 1 lc rgb word(cols,i) noti,\
     for [i=1:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 lc rgb word(cols,i) ti word(hmpk_names,i)
}

if(icomp==3){
plot 1 w l lt -1 noti,\
     for [i=1:words(data_names)] '<paste '.data(word(data_names,i),mesh,snap,fld0,fld0).' '.data(data_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 dt 1 lc i noti,\
     for [i=1:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 lc -1 ti word(hmpk_names,i)
}

}

if(iplot==4){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_residual.eps',name,z)
set output outfile(data_name,z)
}

if(icomp==3){print 'iplot=4 does not work with icomp=3 because the k axis do not align'; print ''; exit}

# Delta^2(k) range
rmin=0.5
rmax=1.5
unset log y
set format y
set yrange [rmin:rmax]
set ylabel 'P_{HM}(k) / P_{OWLS}(k)'

plot NaN w l lw 2 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 2 dt 2 lc -1 ti 'Cross with matter',\
     1 w l lt -1 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld0).' '.hmpk(hmpk_name,z,0,0).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col1 ti 'All matter',\
     '<paste '.data(data_name,mesh,snap,fld1,fld1).' '.hmpk(hmpk_name,z,1,1).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col2 ti 'CDM',\
     '<paste '.data(data_name,mesh,snap,fld2,fld2).' '.hmpk(hmpk_name,z,2,2).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col3 ti 'Gas',\
     '<paste '.data(data_name,mesh,snap,fld3,fld3).' '.hmpk(hmpk_name,z,3,3).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col4 ti 'Stars',\
     '<paste '.data(data_name,mesh,snap,fld6,fld6).' '.hmpk(hmpk_name,z,6,6).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col6 ti 'Pressure',\
     '<paste '.data(data_name,mesh,snap,fld0,fld1).' '.hmpk(hmpk_name,z,0,1).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col2 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld2).' '.hmpk(hmpk_name,z,0,2).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col3 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld3).' '.hmpk(hmpk_name,z,0,3).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col4 noti,\
     '<paste '.data(data_name,mesh,snap,fld0,fld6).' '.hmpk(hmpk_name,z,0,6).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col6 noti

}

if(iplot==5){

# Set fields if none are specified
if(!exists('field1')){field1='all'}
if(!exists('field2')){field2='all'}

# File name for output
if(print==1){
outfile(name,field1,field2)=sprintf('%s_%s_%s_components.eps',name,field1,field2)
set output outfile(data_name,field1,field2)
}

# x axis
set xrange [kmin*1.1:kmax/1.1]
set log x
set xlabel klab

# y axis
plab='{/Symbol D}_{i,j}^2(k)'
pmin=1e-7; pmax=1e3 # Suitable for matter spectra
if(field1 eq 'pressure'){pmin=pmin/1e3; pmax=pmax/1e3}
if(field2 eq 'pressure'){pmin=pmin/1e3; pmax=pmax/1e3}
set yrange [pmin*1.1:pmax/1.1]
set log y
set format y '10^{%T}'
set ylabel plab

# Set field 1
if(field1 eq 'all')     {i1=0}
if(field1 eq 'dm')      {i1=1}
if(field1 eq 'gas')     {i1=2}
if(field1 eq 'stars')   {i1=3}
if(field1 eq 'pressure'){i1=6}

# Set field 2
if(field2 eq 'all')     {i2=0}
if(field2 eq 'dm')      {i2=1}
if(field2 eq 'gas')     {i2=2}
if(field2 eq 'stars')   {i2=3}
if(field2 eq 'pressure'){i2=6}

# Print field information to the screen
print 'Field types:'
print 'all - Matter'
print 'dm - CDM'
print 'gas - Gas'
print 'stars - Stars'
print 'pressure - Pressure'
print 'field 1 *field1*: ', i1, ' ', field1
print 'field 2 *field2*: ', i2, ' ', field2
print ''

unset key

# Title function for plots
title_function(name,field1,field2)=sprintf('Simulation: %s || Power: %s x %s',name,field1,field2)

# y scale info for multiplot
y2=0.9
y1=0.10
my=(y2+y1)/2

# x scale info for multiplot
x1=0.10
x2=0.98
mx=(x2+x1)/2

# Get rid of the generic title
unset title

# Do the actual plotting
set multiplot layout 2,2

do for [i=1:4]{

if(i==1){zz=0.0; ss='snap32'; set xlabel ''; set format x ''; set ylabel plab; set format y '10^{%T}';
set tmargin at screen y2; set lmargin at screen x1; set rmargin at screen mx; set bmargin at screen my;
set label title_function(data_name,field1,field2) at screen 0.37,0.95}

if(i==2){zz=0.5; ss='snap28'; set xlabel ''; set format x ''; set ylabel ''; set format y '';
set tmargin at screen y2; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen my}

if(i==3){zz=1.0; ss='snap26'; set xlabel klab; set format x; set ylabel plab; set format y '10^{%T}';
set tmargin at screen my; set lmargin at screen x1; set rmargin at screen mx; set bmargin at screen y1}

if(i==4){zz=2.0; ss='snap22'; set xlabel klab; set format x; set ylabel ''; set format y '';
set tmargin at screen my; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen y1}

zlab(z)=sprintf('z = %1.1f', z)
set label zlab(zz) at graph 0.05,0.9

plot data(data_name,mesh,ss,field1,field2) u 1:(column(c)-column(s)):5 w e lc 1,\
     hmpk(hmpk_name,zz,i1,i2) u 1:3 w l lc -1 dt 2 lw 3,\
     hmpk(hmpk_name,zz,i1,i2) u 1:4 w l lc -1 dt 3 lw 3,\
     hmpk(hmpk_name,zz,i1,i2) u 1:5 w l lc -1 dt 1 lw 3

unset label

}

unset multiplot

}
