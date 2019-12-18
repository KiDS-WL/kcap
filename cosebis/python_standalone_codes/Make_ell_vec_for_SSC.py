import numpy as np


filename="thps_cov_SSC_KV450_5bins_50_list.dat"
file=open(filename)
list_ell=np.loadtxt(file,comments='#')

ell=list_ell[0:50,9]

filename="input_nonGaussian_ell_vec.ascii"
np.savetxt(filename,ell)

