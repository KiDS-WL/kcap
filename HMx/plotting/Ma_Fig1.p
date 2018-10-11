reset

sun='sun'

file='data/YinZhe_Fig1.dat'
YZ='/Users/Mead/Physics/plots/UPP/UPP_YZM.dat'

set xlabel 'r / Mpc'
set xrange [0:8]

set ylabel 'r^2 P_e(r) / eV cm^{-3} Mpc^2'

set title '10^{15} M_{'.sun.'} / h halo; z = 0.0'

plot file u 1:2 w l lw 3 ti 'UPP from HMx',\
     file u 1:3 w l lw 2 ti 'Pressure profile from HMx',\
     YZ   u 1:2 w l lw 2 ti 'YinZhe'


