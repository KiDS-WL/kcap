unset multiplot
reset

if(!exists('print')){print=0}
if(print==0) {set term aqua dashed size 1000,500}
if(print==1) {set term post enh col dashed dl .5 font ',10'}

print ''

# Plot to make
if(!exists('iplot')){iplot=10}
print 'iplot =  1: OLD: Power spectrum plot'
print 'iplot =  2: OLD: Power spectrum ratio plot'
print 'iplot =  3: Power spectrum suppression plot'
print 'iplot =  4: Power spectrum residual plot'
print 'iplot =  5: Power spectrum components'
print 'iplot =  6: Power spectrum of electron pressure'
print 'iplot =  7: Power spectrum with k^1.5 units'
print 'iplot =  8: PAPER: Combination of iplot=1 and 2'
print 'iplot =  9: PAPER: Response residual'
print 'iplot = 10: Combination of iplot=1 and 2'
print 'iplot = '.iplot.''
print ''

# Simulations to compare against
# 1 - cosmo-OWLS
# 2 - BAHAMAS
# 3 - Generic hydro comparison
if(!exists('icomp')){icomp=3}
if(iplot==8){icomp=3}
print 'icomp = 1: Compare to cosmo-OWLS (NO LONGER SUPPORTED)'
print 'icomp = 2: Compare to BAHAMAS'
print 'icomp = 3: Generic hydro, no comparison simulation'
print 'icomp = '.icomp.''
print ''

if(icomp==1){print 'Twat; icomp=1 does not work'; exit}
# if(icomp==1){sims='cosmo-OWLS'; Om_m=0.272; Om_b=0.0455}
#if(icomp==2){sims='BAHAMAS'}
#if(icomp==3){sims=''}

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
if(iplot==8){z=0.0}
if(z==0.0){snap='snap32'}
if(z==0.5){snap='snap28'}
if(z==1.0){snap='snap26'}
if(z==2.0){snap='snap22'}
snaps="'snap32' 'snap28' 'snap26' 'snap22'"
z_names="'z = 0.0' 'z = 0.5' 'z = 1.0' 'z = 2.0'"
array zs[4]
zs[1]=0.0
zs[2]=0.5
zs[3]=1.0
zs[4]=2.0
print 'Redshift: z = ', z
print ''

# File names - cosmo-OWLS
#if(icomp==1){
#data(sim,type1,type2)=sprintf('/Users/Mead/Physics/cosmo-OWLS/power/N800/%s_%s_%s_power.dat',sim,type1,type2)
#hmpk(sim,i,j)=sprintf('cosmo-OWLS/power_%s_%i%i.dat',sim,i,j)
#print 'Error: cosmo-OWLS comparison needs to be updated'
#exit
#}

if(!exists('mesh')){mesh=1024}
if(iplot==8){mesh=1024}
print 'Mesh size: *mesh*: ', mesh
print ''

# Simulation data files
plot_title_name_z(sim,z)=sprintf('BAHAMAS comarison of %s at z = %1.1f', sim, z)
plot_title_z(z)=sprintf('BAHAMAS comarison at z = %1.1f', z)
simpk(sim,mesh,s,type1,type2)=sprintf('/Users/Mead/Physics/BAHAMAS/power/M%d/%s_nu0_L400N1024_WMAP9_%s_%s_%s_power.dat',mesh,sim,s,type1,type2)
sim_dmonly='DMONLY_2fluid'

# File names - BAHAMAS
if(icomp==2){
hmpk(sim,z,i,j)=sprintf('BAHAMAS/power_%s_z%1.1f_%i%i.dat',sim,z,i,j)
hmpk_dmonly=hmpk('DMONLY',z,0,0)
}

# File names - generic hydro
if(icomp==3){
hmpk(sim,z,i,j)=sprintf('data/power_z%1.1f_%i%i.dat',z,i,j)
hmdm(z)=sprintf('data/power_z%1.1f.dat',z)
hmpk_dmonly=hmdm(z)
}

# Number of columns and pertinent column for simulation power
c=2
s=3
L=5

# Number of columns and pertinent column for halo-model power
d=5
M=5

#cosmo-OWLS simulation names
#if(icomp==1){
#hm_names="'DMONLY' 'REF' 'NOCOOL' 'AGN' 'AGN8p5' 'AGN8p7'"
#owls_names="'DMONLY' 'REF' 'NOCOOL_UVB' 'AGN' 'AGN_Theat_8p5' 'AGN_Theat_8p7'"
#}

# BAHAMAS simulation names
if(icomp==2){hmpk_names="'DMONLY' 'AGN-lo' 'AGN' 'AGN-hi'"}
if(icomp==3){hmpk_names="''"}
sims="'DMONLY_2fluid' 'AGN_7p6' 'AGN_TUNED' 'AGN_8p0'"
sim_names="'DMonly' 'AGN-lo' 'AGN' 'AGN-hi'"

# Set the comparison model
if(!exists('nsim')){nsim=3} # Default to AGN
if(iplot==8){nsim=3} # Default to AGN
hmpk_name=word(hmpk_names,nsim)
print 'Variable: *nsim* '.nsim.''
print 'Halo model power file name: '.hmpk_name.''
sim=word(sims,nsim)
print 'Simuation power file: '.sim.''
print ''

# All different fields for power spectra
fields="'all' 'dm' 'gas' 'stars' '' '' 'epressure'"
field_names="'matter' 'CDM' 'gas' 'stars' '' '' 'electron pressure'"
fld0='all'
fld1='dm'
fld2='gas'
fld3='stars'
fld6='epressure'

print 'Example simulation file: ', simpk(sim,mesh,snap,fld0,fld0)
print 'Example halo-model file: ', hmpk(hmpk_name,z,0,0)
print ''

# Factor to adject the (dimensionful) pressure spectrum
f1=1e2
f2=f1**2
print 'Pressure spectra multiplied by: ', f1
print ''

# Set colours
array cols[6]
cols[1]=1
cols[2]=2
cols[3]=3
cols[4]=4
cols[5]=7
cols[6]=6
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
set ylabel '{/Symbol D}_{uv}^2(k)'
set mytics 10

# Set the overall plot titles
set title plot_title_name_z(sim,z)

## ##

if(iplot==1){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
set output outfile(sim,z)
print 'Outfile: ', outfile(sim,z)
print ''
}

set multiplot layout 1,2

unset title

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     simpk(sim,mesh,snap,fld0,fld0) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col1 noti,\
     simpk(sim,mesh,snap,fld1,fld1) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col2 noti,\
     simpk(sim,mesh,snap,fld2,fld2) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col3 noti,\
     simpk(sim,mesh,snap,fld3,fld3) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col4 noti,\
     simpk(sim,mesh,snap,fld0,fld1) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc col2 noti,\
     simpk(sim,mesh,snap,fld0,fld2) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc col3 noti,\
     simpk(sim,mesh,snap,fld0,fld3) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc col4 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d)) w l lw 3 dt 1 lc col1 ti 'all matter',\
     hmpk(hmpk_name,z,1,1) u 1:(column(d)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(hmpk_name,z,2,2) u 1:(column(d)) w l lw 3 dt 1 lc col3 ti 'gas',\
     hmpk(hmpk_name,z,3,3) u 1:(column(d)) w l lw 3 dt 1 lc col4 ti 'stars',\
     hmpk(hmpk_name,z,0,1) u 1:(column(d)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(hmpk_name,z,0,2) u 1:(column(d)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(hmpk_name,z,0,3) u 1:(column(d)) w l lw 3 dt 2 lc col4 noti

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     simpk(sim,mesh,snap,fld0,fld0) u 1:(column(c)-column(s)):5 w e            pt 7 ps .5 lc col1 noti,\
     simpk(sim,mesh,snap,fld6,fld6) u 1:(f2*(column(c)-column(s))):(f2*$5) w e pt 7 ps .5 lc col6 noti,\
     simpk(sim,mesh,snap,fld0,fld6) u 1:(f1*(column(c)-column(s))):(f1*$5) w e pt 6 ps .5 lc col6 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d))    w l lw 3 dt 1 lc col1 ti 'all matter',\
     hmpk(hmpk_name,z,6,6) u 1:(f2*column(d)) w l lw 3 dt 1 lc col6 ti 'electron pressure',\
     hmpk(hmpk_name,z,0,6) u 1:(f1*column(d)) w l lw 3 dt 2 lc col6 noti

unset multiplot

}

## ##

## ##

if(iplot==2){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_ratio.eps',name,z)
set output outfile(sim,z)
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
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld1,fld1).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col2 noti,\
     '<paste '.simpk(sim,mesh,snap,fld2,fld2).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col3 noti,\
     '<paste '.simpk(sim,mesh,snap,fld3,fld3).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col4 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld1).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col2 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld2).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col3 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld3).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col4 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 ti 'all matter',\
     '<paste '.hmpk(hmpk_name,z,1,1).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     '<paste '.hmpk(hmpk_name,z,2,2).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col3 ti 'gas',\
     '<paste '.hmpk(hmpk_name,z,3,3).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col4 ti 'stars',\
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
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld6).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((f1*column(c))/(column(c+L)-column(s+L)))        w p pt 6 lc col6 noti,\
     '<paste '.simpk(sim,mesh,snap,fld6,fld6).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((f2*column(c))/(column(c+L)-column(s+L)))        w p pt 7 lc col6 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M))    w l lw 3 dt 1 lc col1 ti 'matter',\
     '<paste '.hmpk(hmpk_name,z,0,6).' '.hmpk_dmonly.'' u 1:(f1*column(d)/column(d+M)) w l lw 3 dt 2 lc col6 noti ,\
     '<paste '.hmpk(hmpk_name,z,6,6).' '.hmpk_dmonly.'' u 1:(f2*column(d)/column(d+M)) w l lw 3 dt 1 lc col6 ti 'electron pressure'

unset multiplot

}

## ##

## ##

if(iplot==3){

set title plot_title_z(z)

if(print==1){
outfile='paper/suppression.eps'
set output outfile
}

kmin=0.011
kmax=9.9
set xrange [kmin:kmax]

rmin=0.76
rmax=1.06
unset log y
set format y
set yrange [rmin:rmax]

cols="'black' 'orange' 'orange-red' 'red'"

top=0.98
bot=0.10
lef=0.10
rig=0.98
miy=(top+bot)/2.
mix=(lef+rig)/2.

unset title

set multiplot layout 2,2

do for [i=1:4] {

set key bottom left

if(i==1){snap='snap32'; zlab='z = 0.0'; set format x ''; set xlabel ''; set format y; set ylabel 'P(k) / P_{DMONLY}(k)'}
if(i==2){snap='snap28'; zlab='z = 0.5'; set format x ''; set xlabel ''; set format y ''; set ylabel ''}
if(i==3){snap='snap26'; zlab='z = 1.0'; set format x; set xlabel 'k / h^{-1} Mpc'; set format y; set ylabel 'P(k) / P_{DMONLY}(k)'}
if(i==4){snap='snap22'; zlab='z = 2.0'; set format x; set xlabel 'k / h^{-1} Mpc'; set format y ''; set ylabel ''}

if(i==1){set tmargin at screen top; set bmargin at screen miy; set lmargin at screen lef; set rmargin at screen mix}
if(i==2){set tmargin at screen top; set bmargin at screen miy; set lmargin at screen mix; set rmargin at screen rig}
if(i==3){set tmargin at screen miy; set bmargin at screen bot; set lmargin at screen lef; set rmargin at screen mix}
if(i==4){set tmargin at screen miy; set bmargin at screen bot; set lmargin at screen mix; set rmargin at screen rig}

names=hmpk_names
if(i==2 || i==3 || i==4){names="''''''''"}

set label zlab at graph 0.1,0.9

if(icomp==1 || icomp==2){
plot 1 w l lt -1 noti,\
     for [i=1:words(sims)] '<paste '.simpk(word(sims,i),mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 dt 1 lc rgb word(cols,i) noti,\
     for [i=1:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 lc rgb word(cols,i) ti word(names,i)
}

if(icomp==3){
plot 1 w l lt -1 noti,\
     for [i=1:words(sims)] '<paste '.simpk(word(sims,i),mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 dt 1 lc i noti,\
     for [i=1:words(hmpk_names)] '<paste '.hmpk(word(hmpk_names,i),z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 lc -1 ti word(names,i)
}

unset label

}

unset multiplot

}

## ##

## ##

if(iplot==4){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_residual.eps',name,z)
set output outfile(sim,z)
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
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.hmpk(hmpk_name,z,0,0).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col1 ti 'all matter',\
     '<paste '.simpk(sim,mesh,snap,fld1,fld1).' '.hmpk(hmpk_name,z,1,1).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col2 ti 'CDM',\
     '<paste '.simpk(sim,mesh,snap,fld2,fld2).' '.hmpk(hmpk_name,z,2,2).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col3 ti 'gas',\
     '<paste '.simpk(sim,mesh,snap,fld3,fld3).' '.hmpk(hmpk_name,z,3,3).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col4 ti 'stars',\
     '<paste '.simpk(sim,mesh,snap,fld6,fld6).' '.hmpk(hmpk_name,z,6,6).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 1 lc col6 ti 'electron pressure',\
     '<paste '.simpk(sim,mesh,snap,fld0,fld1).' '.hmpk(hmpk_name,z,0,1).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col2 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld2).' '.hmpk(hmpk_name,z,0,2).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col3 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld3).' '.hmpk(hmpk_name,z,0,3).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col4 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld6).' '.hmpk(hmpk_name,z,0,6).'' u 1:(column(L+d)/(column(c)-column(s))) w l lw 2 dt 2 lc col6 noti

}

## ##

## ##

if(iplot==5){

# Set fields if none are specified
if(!exists('field1')){field1='all'}
if(!exists('field2')){field2='all'}

# File name for output
if(print==1){
outfile(name,field1,field2)=sprintf('%s_%s_%s_components.eps',name,field1,field2)
set output outfile(sim,field1,field2)
}

# x axis
set xrange [kmin*1.1:kmax/1.1]
set log x
set xlabel klab

# y axis
plab='{/Symbol D}_{uv}^2(k)'
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
set label title_function(sim,field1,field2) at screen 0.37,0.95}

if(i==2){zz=0.5; ss='snap28'; set xlabel ''; set format x ''; set ylabel ''; set format y '';
set tmargin at screen y2; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen my}

if(i==3){zz=1.0; ss='snap26'; set xlabel klab; set format x; set ylabel plab; set format y '10^{%T}';
set tmargin at screen my; set lmargin at screen x1; set rmargin at screen mx; set bmargin at screen y1}

if(i==4){zz=2.0; ss='snap22'; set xlabel klab; set format x; set ylabel ''; set format y '';
set tmargin at screen my; set lmargin at screen mx; set rmargin at screen x2; set bmargin at screen y1}

zlab(z)=sprintf('z = %1.1f', z)
set label zlab(zz) at graph 0.05,0.9

plot simpk(sim,mesh,ss,field1,field2) u 1:(column(c)-column(s)):5 w e lc 1,\
     hmpk(hmpk_name,zz,i1,i2) u 1:3 w l lc -1 dt 2 lw 3,\
     hmpk(hmpk_name,zz,i1,i2) u 1:4 w l lc -1 dt 3 lw 3,\
     hmpk(hmpk_name,zz,i1,i2) u 1:5 w l lc -1 dt 1 lw 3

unset label

}

unset multiplot

}

## ##

## ##

if(iplot==6){

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_pressure.eps',name,z)
set output outfile(sim,z)
}

set multiplot layout 2,1

unset xlabel
set format x ''

set key top left

plot simpk(sim,mesh,snap,fld0,fld0) u 1:(column(c)-column(s)):5            w e pt 7 lc col1 noti,\
     simpk(sim,mesh,snap,fld0,fld6) u 1:(f1*(column(c)-column(s))):(f1*$5) w e pt 6 lc col5 noti,\
     simpk(sim,mesh,snap,fld6,fld6) u 1:(f2*(column(c)-column(s))):(f2*$5) w e pt 7 lc col6 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d))      w l lw 3 dt 1 lc col1 ti 'matter-matter',\
     hmpk(hmpk_name,z,0,0) u 1:(column(d-2))    w l lw 3 dt 2 lc col1 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d-1))    w l lw 3 dt 3 lc col1 noti,\
     hmpk(hmpk_name,z,0,6) u 1:(f1*column(d))   w l lw 3 dt 1 lc col5 ti 'matter-electron pressure',\
     hmpk(hmpk_name,z,0,6) u 1:(f1*column(d-2)) w l lw 3 dt 2 lc col5 noti,\
     hmpk(hmpk_name,z,0,6) u 1:(f1*column(d-1)) w l lw 3 dt 3 lc col5 noti,\
     hmpk(hmpk_name,z,6,6) u 1:(f2*column(d))   w l lw 3 dt 1 lc col6 ti 'electron pressure-electron ressure',\
     hmpk(hmpk_name,z,6,6) u 1:(f2*column(d-2)) w l lw 3 dt 2 lc col6 noti,\
     hmpk(hmpk_name,z,6,6) u 1:(f2*column(d-1)) w l lw 3 dt 3 lc col6 noti

unset title

set xlabel klab
set format x

# Set the y axis for P(k) / P_{DMONLY}(k)
rmin=1e-5
rmax=2
set yrange [rmin:rmax]
set ylabel 'P(k) / P_{DMONLY}(k)'

plot 1 w l lt -1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L)))    w p pt 7 lc col1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld6).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:(f1*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 lc col5 noti,\
     '<paste '.simpk(sim,mesh,snap,fld6,fld6).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:(f2*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 lc col6 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M))    w l lw 3 dt 1 lc col1 noti '{/Symbol d}{/Symbol d}',\
     '<paste '.hmpk(hmpk_name,z,0,6).' '.hmpk_dmonly.'' u 1:(f1*column(d)/column(d+M)) w l lw 3 dt 1 lc col5 noti '{/Symbol d}p',\
     '<paste '.hmpk(hmpk_name,z,6,6).' '.hmpk_dmonly.'' u 1:(f2*column(d)/column(d+M)) w l lw 3 dt 1 lc col6 noti 'pp'

unset multiplot

}

## ##

## ##

if(iplot==7){

pow=1.5

# Delta^2(k)/k^1.5 range
pmin=1e-4
pmax=1e2
set log y
set yrange [pmin:pmax]
set format y '10^{%T}'
set ylabel '{/Symbol D}_{uv}^2(k) / [k / h^{-1} Mpc]^{1.5}'
set mytics 10

if(print==1){
outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
set output outfile(sim,z)
}

set key bottom right

plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     simpk(sim,mesh,snap,fld0,fld0) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col1 noti,\
     simpk(sim,mesh,snap,fld1,fld1) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col2 noti,\
     simpk(sim,mesh,snap,fld2,fld2) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col3 noti,\
     simpk(sim,mesh,snap,fld3,fld3) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 7 lc col4 noti,\
     simpk(sim,mesh,snap,fld0,fld1) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 6 lc col2 noti,\
     simpk(sim,mesh,snap,fld0,fld2) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 6 lc col3 noti,\
     simpk(sim,mesh,snap,fld0,fld3) u 1:((column(c)-column(s))/(column(1)**pow)):($5/(column(1)**pow)) w e pt 6 lc col4 noti,\
     hmpk(hmpk_name,z,0,0) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col1 ti 'all matter',\
     hmpk(hmpk_name,z,1,1) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(hmpk_name,z,2,2) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col3 ti 'gas',\
     hmpk(hmpk_name,z,3,3) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 1 lc col4 ti 'stars',\
     hmpk(hmpk_name,z,0,1) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(hmpk_name,z,0,2) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(hmpk_name,z,0,3) u 1:((column(d)/column(1)**pow)) w l lw 3 dt 2 lc col4 noti

}

## ##

## PAPER FIGURE: big ##

if(iplot==8 || iplot==10){

if(print==0){set term aqua dashed size 1200,800}

if(print==1){
#outfile(name,z)=sprintf('%s_z%1.1f_power.eps',name,z)
outfile='paper/hydro.eps'
set output outfile
print 'Outfile: ', outfile
print ''
}

set multiplot layout 2,2

set label ''.word(sim_names,nsim).'; z = '.sprintf('%1.1f', z).'' at graph 0.05,0.9

set xlabel ''
set format x ''

set key bottom right

unset title

# Top left - matter spectra
plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     simpk(sim,mesh,snap,fld0,fld0) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col1 noti,\
     simpk(sim,mesh,snap,fld1,fld1) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col2 noti,\
     simpk(sim,mesh,snap,fld2,fld2) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col3 noti,\
     simpk(sim,mesh,snap,fld3,fld3) u 1:(column(c)-column(s)):5 w e pt 7 ps .5 lc col4 noti,\
     simpk(sim,mesh,snap,fld0,fld1) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc col2 noti,\
     simpk(sim,mesh,snap,fld0,fld2) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc col3 noti,\
     simpk(sim,mesh,snap,fld0,fld3) u 1:(column(c)-column(s)):5 w e pt 6 ps .5 lc col4 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d)) w l lw 3 dt 1 lc col1 ti 'all matter',\
     hmpk(hmpk_name,z,1,1) u 1:(column(d)) w l lw 3 dt 1 lc col2 ti 'CDM',\
     hmpk(hmpk_name,z,2,2) u 1:(column(d)) w l lw 3 dt 1 lc col3 ti 'gas',\
     hmpk(hmpk_name,z,3,3) u 1:(column(d)) w l lw 3 dt 1 lc col4 ti 'stars',\
     hmpk(hmpk_name,z,0,1) u 1:(column(d)) w l lw 3 dt 2 lc col2 noti,\
     hmpk(hmpk_name,z,0,2) u 1:(column(d)) w l lw 3 dt 2 lc col3 noti,\
     hmpk(hmpk_name,z,0,3) u 1:(column(d)) w l lw 3 dt 2 lc col4 noti

unset label

# Top right - pressure spectra
plot NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     simpk(sim,mesh,snap,fld0,fld0) u 1:(column(c)-column(s)):5            w e pt 7 ps .5 lc col1 noti,\
     simpk(sim,mesh,snap,fld6,fld6) u 1:(f2*(column(c)-column(s))):(f2*$5) w e pt 7 ps .5 lc col6 noti,\
     simpk(sim,mesh,snap,fld0,fld6) u 1:(f1*(column(c)-column(s))):(f1*$5) w e pt 6 ps .5 lc col6 noti,\
     hmpk(hmpk_name,z,0,0) u 1:(column(d))    w l lw 3 dt 1 lc col1 ti 'all matter',\
     hmpk(hmpk_name,z,6,6) u 1:(f2*column(d)) w l lw 3 dt 1 lc col6 ti 'electron pressure',\
     hmpk(hmpk_name,z,0,6) u 1:(f1*column(d)) w l lw 3 dt 2 lc col6 noti

set xlabel 'k / h Mpc^{-1}'
set format x

rmin=3e-5
rmax=2
set ylabel 'P_{uv}(k) / P_{no-hydro}(k)'
set yrange [rmin:rmax]

# Bottom left - matter response
plot 1 w l lt -1 noti,\
     Om_b/Om_m      w l lc -1 dt 2 noti,\
     Om_c/Om_m      w l lc -1 dt 2 noti,\
     (Om_b/Om_m)**2 w l lc -1 dt 2 noti,\
     (Om_c/Om_m)**2 w l lc -1 dt 2 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc col1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld1,fld1).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc col2 noti,\
     '<paste '.simpk(sim,mesh,snap,fld2,fld2).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc col3 noti,\
     '<paste '.simpk(sim,mesh,snap,fld3,fld3).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc col4 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 ps .5 lc col1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld1).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 ps .5 lc col2 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld2).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 ps .5 lc col3 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld3).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 ps .5 lc col4 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col1 noti,\
     '<paste '.hmpk(hmpk_name,z,1,1).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col2 noti,\
     '<paste '.hmpk(hmpk_name,z,2,2).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col3 noti,\
     '<paste '.hmpk(hmpk_name,z,3,3).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 1 lc col4 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col1 noti,\
     '<paste '.hmpk(hmpk_name,z,0,1).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col2 noti,\
     '<paste '.hmpk(hmpk_name,z,0,2).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col3 noti,\
     '<paste '.hmpk(hmpk_name,z,0,3).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M)) w l lw 3 dt 2 lc col4 noti

# Bottom right - pressure response
plot 1 w l lt -1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld0).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:((column(c)-column(s))/(column(c+L)-column(s+L)))    w p pt 7 ps .5 lc col1 noti,\
     '<paste '.simpk(sim,mesh,snap,fld6,fld6).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:(f2*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 7 ps .5 lc col6 noti,\
     '<paste '.simpk(sim,mesh,snap,fld0,fld6).' '.simpk(sim_dmonly,mesh,snap,fld0,fld0).'' u 1:(f1*(column(c)-column(s))/(column(c+L)-column(s+L))) w p pt 6 ps .5 lc col6 noti,\
     '<paste '.hmpk(hmpk_name,z,0,0).' '.hmpk_dmonly.'' u 1:(column(d)/column(d+M))    w l lw 3 dt 1 lc col1 noti,\
     '<paste '.hmpk(hmpk_name,z,6,6).' '.hmpk_dmonly.'' u 1:(f2*column(d)/column(d+M)) w l lw 3 dt 1 lc col6 noti,\
     '<paste '.hmpk(hmpk_name,z,0,6).' '.hmpk_dmonly.'' u 1:(f1*column(d)/column(d+M)) w l lw 3 dt 2 lc col6 noti

unset multiplot

}

## ##

## PAPER: Response ratio ##

if(iplot==9){

# Set to do BAHAMAS comparison only
icomp=2
print 'Note that icomp=2 has been automatically set here'
print ''

# Terminal commands
if(print==0){set term aqua dashed size 1200,800}
if(print==1){set term post enh col dashed; set output 'paper/response_ratio.eps'}

# Label position
labx=0.25
laby=0.9

# Plot limits in x
lef=0.08
rig=0.98
nx=3.
dx=(rig-lef)/nx

# k axis
kmin=0.013
kmax=9
set xrange [kmin:kmax]
set log x

# Plot limits in y
top=0.99
bot=0.08
ny=4.
dy=(top-bot)/ny

# ratio axis
rmin=0.5
rmax=1.5
unset log y
set yrange [rmin:rmax]

# No title for this plot
unset title

set key top left

combi(isim,iz,f1,f2)=sprintf('<paste '.hmpk(word(hmpk_names,isim),zs[iz],f1,f2).' '.hmpk('DMONLY',zs[iz],0,0).' '.data(word(sims,isim),mesh,word(snaps,iz),word(fields,f1+1),word(fields,f2+1)).' '.data('DMONLY_2fluid',mesh,word(snaps,iz),fld0,fld0).'',isim,iz,f1,f2)

isim=3
iz=2

set multiplot layout 4,3

do for [isim=2:4]{

if(isim==2){set lmargin at screen lef+0*dx; set rmargin at screen lef+1*dx}
if(isim==3){set lmargin at screen lef+1*dx; set rmargin at screen lef+2*dx}
if(isim==4){set lmargin at screen lef+2*dx; set rmargin at screen lef+3*dx}

if(isim==2){set format y; set ylabel 'R_{HM} / R_{sim}'}
if(isim==3 || isim==4){set format y ''; set ylabel ''}

do for [iz=1:4]{

if(iz==1){set tmargin at screen top-0*dy; set bmargin at screen top-1*dy}
if(iz==2){set tmargin at screen top-1*dy; set bmargin at screen top-2*dy}
if(iz==3){set tmargin at screen top-2*dy; set bmargin at screen top-3*dy}
if(iz==4){set tmargin at screen top-3*dy; set bmargin at screen top-4*dy}

if(iz==1 || iz==2 || iz==3){set xlabel ''; set format x ''; unset key}
if(iz==4){set xlabel 'k / h^{-1} Mpc'; set format x}

set label ''.word(hmpk_names,isim).'; '.word(z_names,iz).'' at graph labx,laby

if(iz==1 && isim==2){set key; unset label}

plot 1 w l lt -1 noti,\
     NaN w l lw 3 dt 1 lc -1 ti 'Autospectra',\
     NaN w l lw 3 dt 2 lc -1 ti 'Cross with matter',\
     combi(isim,iz,0,0) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[1] ti word(field_names,1+0),\
     combi(isim,iz,1,1) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[2] ti word(field_names,1+1),\
     combi(isim,iz,2,2) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[3] ti word(field_names,1+2),\
     combi(isim,iz,3,3) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[4] ti word(field_names,1+3),\
     combi(isim,iz,6,6) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 1 lc cols[6] ti word(field_names,1+6),\
     combi(isim,iz,0,1) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[2] noti 'CDM',\
     combi(isim,iz,0,2) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[3] noti 'gas',\
     combi(isim,iz,0,3) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[4] noti 'stars',\
     combi(isim,iz,0,6) u 1:((column(d)/column(d+M))/((column(2*M+c)-column(2*M+s))/(column(2*M+c+L)-column(2*M+s+L)))) w l lw 3 dt 2 lc cols[6] noti 'electron_pressure'

unset label

}

}

unset multiplot

}
