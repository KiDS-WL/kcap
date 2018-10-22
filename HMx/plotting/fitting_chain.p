reset

chain='fitting/test_chain.dat'

np=12

plot for [i=1:np] chain u (column(0)):(column(i+1)) w p noti
