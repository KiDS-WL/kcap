import numpy as np


def cutCov(Cov_in,nPairs,nbins_in,nbins_min,nbins_max):
	nbins=nbins_max-nbins_min
	CovCut=np.zeros((nbins*nPairs,nbins*nPairs))
	for i1 in range(nPairs):
		for i2 in range(nPairs):
			c1_min=i1*nbins
			c2_min=i2*nbins
			c1_max=i1*nbins+nbins
			c2_max=i2*nbins+nbins
	# 
			i1_min=i1*nbins_in+nbins_min
			i2_min=i2*nbins_in+nbins_min
			i1_max=i1*nbins_in+nbins_max
			i2_max=i2*nbins_in+nbins_max
			CovCut[c1_min:c1_max,c2_min:c2_max]=Cov_in[i1_min:i1_max,i2_min:i2_max]
	return CovCut



def cutData(Data_in,nPairs,nbins_in,nbins_min,nbins_max):
	nbins=nbins_max-nbins_min
	Data_cut=np.zeros((nbins*nPairs))
	for i1 in range(nPairs):
		c1_min=i1*nbins
		c1_max=i1*nbins+nbins
	# 
		i1_min=i1*nbins_in+nbins_min
		i1_max=i1*nbins_in+nbins_max
		Data_cut[c1_min:c1_max]=Data_in[i1_min:i1_max]
	return Data_cut

def cutCov2PCFs(Cov_in,nPairs,nTheta,nTheta_plus_min,nTheta_plus_max,nTheta_minus_min,nTheta_minus_max):
	nTheta_plus=nTheta_plus_max-nTheta_plus_min
	nTheta_minus=nTheta_minus_max-nTheta_minus_min
	CovCut=np.zeros(((nTheta_plus+nTheta_minus)*nPairs,(nTheta_plus+nTheta_minus)*nPairs))
	for i1 in range(nPairs):
		for i2 in range(nPairs):
			# first for ++
			c1_min=i1*nTheta_plus
			c2_min=i2*nTheta_plus
			c1_max=i1*nTheta_plus+nTheta_plus
			c2_max=i2*nTheta_plus+nTheta_plus
	# 
			i1_min=i1*nTheta+nTheta_plus_min
			i2_min=i2*nTheta+nTheta_plus_min
			i1_max=i1*nTheta+nTheta_plus_max
			i2_max=i2*nTheta+nTheta_plus_max
			CovCut[c1_min:c1_max,c2_min:c2_max]=Cov_in[i1_min:i1_max,i2_min:i2_max]
			# now for --
			c1_min=i1*nTheta_minus+nPairs*nTheta_plus
			c2_min=i2*nTheta_minus+nPairs*nTheta_plus
			c1_max=i1*nTheta_minus+nTheta_minus+nPairs*nTheta_plus
			c2_max=i2*nTheta_minus+nTheta_minus+nPairs*nTheta_plus
			# print(c1_min,c1_max,c2_min,c2_max)
	# 
			i1_min=i1*nTheta+nTheta_minus_min+nPairs*nTheta
			i2_min=i2*nTheta+nTheta_minus_min+nPairs*nTheta
			i1_max=i1*nTheta+nTheta_minus_max+nPairs*nTheta
			i2_max=i2*nTheta+nTheta_minus_max+nPairs*nTheta
			# print(i1_min,i1_max,i2_min,i2_max)
			CovCut[c1_min:c1_max,c2_min:c2_max]=Cov_in[i1_min:i1_max,i2_min:i2_max]
			# now for -+
			c1_min=i1*nTheta_minus+nPairs*nTheta_plus
			c1_max=i1*nTheta_minus+nTheta_minus+nPairs*nTheta_plus
			c2_min=i2*nTheta_plus
			c2_max=i2*nTheta_plus+nTheta_plus
	# 
			i1_min=i1*nTheta+nTheta_minus_min+nPairs*nTheta
			i2_min=i2*nTheta+nTheta_plus_min
			i1_max=i1*nTheta+nTheta_minus_max+nPairs*nTheta
			i2_max=i2*nTheta+nTheta_plus_max
			CovCut[c1_min:c1_max,c2_min:c2_max]=Cov_in[i1_min:i1_max,i2_min:i2_max]
			# now for +-
			c1_min=i1*nTheta_plus
			c2_min=i2*nTheta_minus+nPairs*nTheta_plus
			c1_max=i1*nTheta_plus+nTheta_plus
			c2_max=i2*nTheta_minus+nTheta_minus+nPairs*nTheta_plus
	# 
			i1_min=i1*nTheta+nTheta_plus_min
			i2_min=i2*nTheta+nTheta_minus_min+nPairs*nTheta
			i1_max=i1*nTheta+nTheta_plus_max
			i2_max=i2*nTheta+nTheta_minus_max+nPairs*nTheta
			CovCut[c1_min:c1_max,c2_min:c2_max]=Cov_in[i1_min:i1_max,i2_min:i2_max]
	return CovCut


# under construction
def cutCovMultiple(Cov_in,nPairs,nSections,nbins_in,nbins_min,nbins_max):
	nbins=nbins_min-nbins_max
	# check if an array input is given 
	if hasattr(nbins, "__len__"):
		nbins_tot=nbins.sum()
	else:
		nbins_tot=nbins
	# 
	CovCut=np.zeros((nbins_tot*nPairs,nbins_tot*nPairs))
	nbins_sum=0
	# first do the auto correlations
	for ibins in range(nbins_tot):
		start=nbins_sum*nPairs
		end  =(nbins_sum+nbins_in[ibins])*nPairs
		Cov_ibins_in=Cov_in[start:end,start:end]
		CovCut[start:end,start:end]=cutCov(Cov_ibins_in,nPairs[ibins],nbins_in[ibins],nbins_min[ibins],nbins_max[ibins])
		nbins_sum+=nbins_in[ibins]

	# now the cross-correlations, not done yet
	for i1 in range(nPairs):
		for i2 in range(nPairs):
			# now for 1 and 2
			c1_min=i1*nTheta_minus+nPairs*nTheta_plus
			c1_max=i1*nTheta_minus+nTheta_minus+nPairs*nTheta_plus
			c2_min=i2*nTheta_plus
			c2_max=i2*nTheta_plus+nTheta_plus
	# 
			i1_min=i1*nTheta+nTheta_minus_min+nPairs*nTheta
			i2_min=i2*nTheta+nTheta_plus_min
			i1_max=i1*nTheta+nTheta_minus_max+nPairs*nTheta
			i2_max=i2*nTheta+nTheta_plus_max
			CovCut[c1_min:c1_max,c2_min:c2_max]=Cov_in[i1_min:i1_max,i2_min:i2_max]
			# now for 2 and 1
			c1_min=i1*nTheta_plus
			c2_min=i2*nTheta_minus+nPairs*nTheta_plus
			c1_max=i1*nTheta_plus+nTheta_plus
			c2_max=i2*nTheta_minus+nTheta_minus+nPairs*nTheta_plus
	# 
			i1_min=i1*nTheta+nTheta_plus_min
			i2_min=i2*nTheta+nTheta_minus_min+nPairs*nTheta
			i1_max=i1*nTheta+nTheta_plus_max
			i2_max=i2*nTheta+nTheta_minus_max+nPairs*nTheta
			CovCut[c1_min:c1_max,c2_min:c2_max]=Cov_in[i1_min:i1_max,i2_min:i2_max]

	return CovCut

# cut redshift pairs in the cov
# under construction
def cutCovRedshift(Cov_in,nPairs,nSections,nbins_in,pairs_to_remove):
	# remove the pairs_to_remove from a cross covariance with nSections and nPairs, 
	# each section has nbins_in x_bins.
	if hasattr(pairs_to_remove, "__len__"):
		nPairs_remove=len(pairs_to_remove)
	else:
		nPairs_remove=1
	
	CovCut=np.zeros((nbins_in*(nPairs-nPairs_remove),nbins_tot*(nPairs-nPairs_remove)))
	#This needs to be written

	return cutCov


