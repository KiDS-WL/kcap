reset

if(print==0){set term aqua dashed}
if(print==1){set term post enh col; set output 'HMx_params.eps'}

data(name)=sprintf('data/HMx_params_%s.dat',name)

names='AGN-lo AGN AGN-hi'

set multiplot layout 3,2

set xlabel 'z'

do for [j=1:6]{

if(j==1){ylab='{/Symbol a}'}
if(j==2){ylab='{/Symbol e}'}
if(j==3){ylab='{/Symbol G}'}
if(j==4){ylab='M_0'}
if(j==5){ylab='A_*'}
if(j==6){ylab='T_w'}

set ylabel ylab

plot for [i=1:3] NaN w l lc i dt 1 lw 3 ti word(names,i),\
     for [i=1:3] data(word(names,i)) u 1:(column(1+j)) w l lc i dt 1 lw 3 noti

}

unset multiplot
