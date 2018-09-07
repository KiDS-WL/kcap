reset

if(print==0) set term aqua dashed
if(print==1) set term post enh col; set output 'CFHTLenS.eps'

types="'linear' '2halo' '1halo' 'full'"
xi(type)=sprintf('data/xi_%s.dat',type)

data='/Users/Mead/Physics/CFHTLenS/corr.txt'

xip_col='red'
xim_col='orange'

xip_lab='{/Symbol x}_+'
xim_lab='{/Symbol x}_-'

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

set log y
set ylabel '{/Symbol x}({/Symbol q})'
set format y '10^{%T}'
set yrange [1e-8:1e-4]
set mytics 10

plot small w l lw 3 lc rgb xip_col dt 1 ti xip_lab,\
     small w l lw 3 lc rgb xim_col dt 1 ti xim_lab,\
     small w l lw 3 lc -1 dt 1 ti 'Full',\
     small w l lw 3 lc -1 dt 2 ti 'One halo',\
     small w l lw 3 lc -1 dt 3 ti 'Two halo',\
     small w l lw 3 lc -1 dt 4 ti 'Linear',\
     for [i=1:words(types)] xi(word(types,i)) u 1:2 w l lw 3 lc rgb 'red'    dt 5-i noti,\
     for [i=1:words(types)] xi(word(types,i)) u 1:4 w l lw 3 lc rgb 'orange' dt 5-i noti,\
     data u ($1/arcmin):2 w p ps 1.5 pt 7 lc rgb 'black' ti 'CFHTLenS {/Symbol x}_+',\
     data u ($1/arcmin):2:(sqrt($3)) w errorbars ps 0 lt 1 lc rgb 'black' noti,\
     data u (offset*$1/arcmin):4 w p ps 1.5 pt 7 lc rgb 'dark-grey' ti 'CFHTLenS {/Symbol x}_-',\
     data u (offset*$1/arcmin):4:(sqrt($5)) w errorbars ps 0 lt 1 lc rgb 'dark-grey' noti
