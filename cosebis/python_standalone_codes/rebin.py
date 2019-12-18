import numpy as np

# output is x_log_binned, x_weighted_binned, signal_weighted_binned
def rebin(x,signal,weight,x_min,x_max,nbins):
	# print('rebinning now')
	binned_output=np.zeros((nbins,3))
	for ibins in range(nbins):
		x_binned=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+0.5))
		upperEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+1.0))
		lowerEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins))
		good=((x<upperEdge) & (x>lowerEdge))
		# print(x_binned)
		if(good.any()):
			weight_sum=weight[good].sum()
			x_binned_weighted=(x[good]*weight[good]).sum()/weight_sum
			binned_output[ibins,0]=x_binned
			binned_output[ibins,1]=x_binned_weighted
			binned_output[ibins,2]=(signal[good]*weight[good]).sum()/weight_sum
			# print(ibins,weight_sum,len(weight[good]))
		else:
			print("WARNING: not enough bins to rebin to "+str(nbins)+" log bins")
	return binned_output

# nbins_in=300
# x_min=0.25
# x_max=400.0
# nbins=8
# ibins=np.linspace(0,nbins_in-1,nbins_in)
# x=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins_in)*(ibins+0.5))
# weight=np.ones(len(x))

# upperEdge_in=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins_in)*(ibins+1.0))
# lowerEdge_in=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins_in)*(ibins))
# Deltax_in=np.log(upperEdge_in/lowerEdge_in)

# ibins=np.linspace(0,nbins-1,nbins)
# x_binned=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+0.5))
# # Deltax=np.ones(nbins_in)
# for ibins in range(nbins):
# 	upperEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+1.0))
# 	lowerEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins))
# 	Deltax=np.log(upperEdge/lowerEdge)
# 	good=((x<upperEdge) & (x>lowerEdge))
# 	x_binned_log=np.exp(np.log(x[good]).sum()/len(x[good]))
# 	x_binned_weighted=(x[good]*weight[good]).sum()/weight[good].sum()
# 	# x_binned=np.exp(np.log(x[good]).sum()/float(len(x[good])))
# 	print(x_binned[ibins],x_binned_log,x_binned_weighted)

