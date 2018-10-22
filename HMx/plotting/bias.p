unset multiplot
reset

if(!exists('print')){print=0}
if(print==0){set term aqua dashed}
if(print==1){set term post enh col sol; set output 'bias.eps'}

kmin=1e-3
kmax=1e2
set log x
set xrange [kmin:kmax]

pmin=1e-7
pmax=1e4
set log y
set yrange [pmin:pmax]
#set yrange [*:*]
set ylabel '{/Symbol D}_{ij}^2(k)'
set format y '10^{%T}'

#File to compare against
if(!exists('com')){com='power_HMcode_standard.txt'}

# Number of redshifts
nz=16

power(a,b,type)=sprintf('data/power_%s%s_%s.dat',a,b,type)

print ''
print 'Comparing power_hm.dat to: ', com
print 'Set variable *com* as comparison file'
print ''

set key top left

col_red(i)=sprintf("#%1x%1x0000",i-1,i-1)
col_green(i)=sprintf("#00%1x%1x00",i-1,i-1)

set lmargin 10
set rmargin 2

set multiplot layout 2,1

set xlabel ''
set format x ''

plot for[i=1:1]  power('m','m','hm') u 1:(column(i+1)) w l lw 2 dt 1 lc -1 ti 'Matter-Matter',\
     for[i=1:1]  power('m','f','hm') u 1:(column(i+1)) w l lw 2 dt 1 lc 1  ti 'Field-Matter',\
     for[i=1:1]  power('f','f','hm') u 1:(column(i+1)) w l lw 2 dt 1 lc 2  ti 'Field-Field'

set xlabel 'k / (h Mpc^{-1})'
set format x

unset log y
set yrange [0:2]
set format y
set ylabel 'b_f(k)'

plot 1 w l ls -1 noti,\
     for[i=1:1] '<paste '.power('m','f','hm').' '.power('m','m','hm').'' u 1:(column(i+1)/column(i+2+nz))       w l lw 2 dt 1 lc 1 noti,\
     for[i=1:1] '<paste '.power('f','f','hm').' '.power('m','m','hm').'' u 1:(sqrt(column(i+1)/column(i+2+nz))) w l lw 2 dt 1 lc 2 noti

unset multiplot



