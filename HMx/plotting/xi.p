reset

if(print==0) set term aqua
if(print==1) set term post enh col; set output 'xi.eps'

types="'linear' '2halo' '1halo' 'full'"
xi(type)=sprintf('data/xi_%s.dat',type)

small=1e-18
arcmin=60.
offset=0.95

thmin=0.01
thmax=10.

set log x
set xlabel '{/Symbol q} [degrees]'
set xrange [thmin:thmax]
set xtics nomirror

#Set x2 axis for arcminutes
set log x2
set format x2
set x2range [thmin*arcmin:thmax*arcmin]
set x2label '{/Symbol q} [arcminutes]'
set x2tics
set mx2tics 10

set log y
set ylabel '{/Symbol x}({/Symbol q})'
set format y '10^{%T}'
#set yrange [1e-8:1e-4]
set yrange [1e-12:1e-8]

plot small w l lw 3 dt 1 lc 1  ti '{/Symbol x} with J_0',\
     small w l lw 3 dt 1 lc 2  ti '{/Symbol x} with J_2',\
     small w l lw 3 dt 1 lc 3  ti '{/Symbol x} with J_4',\
     small w l lw 3 dt 1 lc -1 ti 'Full',\
     small w l lw 3 dt 2 lc -1 ti 'One halo',\
     small w l lw 3 dt 3 lc -1 ti 'Two halo',\
     small w l lw 3 dt 4 lc -1 ti 'Linear',\
for [j=1:3] for [i=1:words(types)] xi(word(types,i)) u 1:(column(j+1)) w l lw 3 lc j    dt 5-i noti
