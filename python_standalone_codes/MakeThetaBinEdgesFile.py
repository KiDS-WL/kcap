import numpy as np
import matplotlib.pylab as plt

nbins=9
x_min=0.5
x_max=300.0

ix=np.linspace(0,nbins-1,nbins)
# bin mid point for x
# x_mid=np.exp(np.log(x_min)+(np.log(x_max)-np.log10(x_min))/(nbins_broad)*(ix+0.5));
x=np.exp(np.log(x_min)+(np.log(x_max)-np.log(x_min))/(nbins)*(ix))
x_min_edge=np.exp(np.log(x_min)+(np.log(x_max)-np.log(x_min))/(nbins)*(ix))
x_max_edge=np.exp(np.log(x_min)+(np.log(x_max)-np.log(x_min))/(nbins)*(ix+1.0))


x_edges=np.transpose(np.vstack((x_min_edge,x_max_edge)))

filename="../example_files/xi/inputs/theta_bin_edges_file.ascii"
np.savetxt(filename,x_edges,comments="theta bin edges")