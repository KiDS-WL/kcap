reset

kmin=0.02
kmax=10
set log x
set xlabel 'k / h Mpc^{-1}'
set xrange [kmin:kmax]

dr=0.2
rmin=0.9
rmax=1.16
set yrange[rmin:rmax]
set ylabel 'P(k) / P_{{/Symbol L}CDM}(k)'

power(name)=sprintf('data/power_%s.dat',name)
ratio(name)=sprintf('/Users/Mead/Physics/fixed_sigma/power/L200/z0/%s_ratio.dat',name)

names="'LCDM' 'OCDM' 'w-0.7' 'w-1.3' 'wa0.5' 'wa-0.5' 'w0-0.7wa-1.5' 'w0-1.3wa0.5'"
titles="'{/Symbol L}CDM' 'Open' 'w = -0.7' 'w = -1.3' 'w_a = 0.5' 'w_a = -0.5' 'w = -0.7; w_a = -1.5' 'w = -1.3; w_a = 0.5'"

set key top left

plot 1 w l lt -1 noti,\
     for [i=2:words(names)] '<paste '.power(word(names,i)).' '.power('LCDM').'' u 1:($5/$10) w l lc i lw 2 ti word(titles,i),\
     for [i=2:words(names)] ratio(word(names,i)) u 1:2:3 w errorbars pt 7 ps .5 dt 1 lc i noti
