# This module computes the average luminosity scaling of the intrinsic alignment 
# following the formalism presented in Joachimi et al. 2011b
# 
#							A(L) = (L/L_0)^beta
#
# If the luminosity distribution of the sample is not narrow, we need to average
# the contribution from L^beta of the entire sample, i.e.
#
#			<(L/L_0)^beta> = (1/L_0^beta) \int L^beta p(L) dL
#
# where p(L) is the pdf of L.


# -------------------------------------------------------------------------------- #
# IMPORTANT: here we assume luminosities to be in units of L_sun/h2
# -------------------------------------------------------------------------------- #


import sys
import re
from os import listdir
import numpy as np
import clf_lib as lf
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import simps
import matplotlib.pyplot as plt

import time



# library
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
    
    

def load_data(file_name):
	z_data, loglum_data = np.loadtxt(file_name,unpack=True, dtype=np.float)
	return z_data, loglum_data
	
def create_z_bins(zmin, zmax, nz):
	z_histo = np.linspace(zmin, zmax, nz+1)
	dz = (zmax-zmin)/nz
	z_cen = z_histo[0:-1] + 0.5*dz
	return z_cen, z_histo, dz
	
def luminosity_pdf(xlum, loglum_data, nlbins):
    lum = 10.**loglum_data
    lum_hist, lum_bins = np.histogram(lum, bins=nlbins)
    dbin = (lum_bins[-1]-lum_bins[0])/(nlbins)
    bincen = lum_bins[0:-1]+0.5*dbin
    norm_pdf = simps(lum_hist, bincen)
    pdf = lum_hist/norm_pdf
    print simps(pdf,bincen)
    interp = interp1d(bincen, pdf, fill_value=(0,0), bounds_error=False)
    print 'extrapolation worked fine'
    lum_pdf = interp(xlum)   
    plt.loglog(xlum, lum_pdf, 'r-')
    plt.loglog(bincen, pdf, 'b--')

    plt.show()
    return lum_pdf

def configure():
	# ------------------------------------ input section ------------------------------------------ #
	directory_name = "/Users/mcfortuna/IA/data/MICE/mice_1000deg2_fluxlim_r24/6_bins/red_cen_lum/" 
	outfile = "/Users/mcfortuna/IA/data/MICE/mice_1000deg2_fluxlim_r24/6_bins/red_cen_lum_pdf.txt"
	# --------------------------------------------------------------------------------------------- #
	ldir = listdir(directory_name)
	binfile = [i for i in ldir if i.startswith('lum_')]
	binfile.sort(key=natural_keys)
	nfiles = len(binfile)
	print('nfiles = ', nfiles)	
	print(binfile)
	zbin = {}	
	loglum = {}
	for i in range(nfiles):
		zbin["z{0}".format(i+1)], loglum["z{0}".format(i+1)] = load_data(directory_name+binfile[i])		
	return zbin, loglum, outfile



def main():	
	zbin, loglum, outfile = configure()	
	nz = len(zbin)		
	nlum = 2000
	lum_pdf = np.empty([nz,nlum])
	lum_xarray = np.logspace(7., 12., nlum)
	
	for jz in range(0,nz):
		z_j = 'z%s' %str(jz+1)
		lum_pdf[jz] = luminosity_pdf(lum_xarray, loglum[z_j], nlum)
		
	# concatenate the bins for the txt file
	lum_pdf_binz = np.c_[lum_xarray, lum_pdf[0]]
	
	for jz in range(1,nz):
		lum_pdf_binz = np.c_[lum_pdf_binz, lum_pdf[jz]]

	np.savetxt(outfile, lum_pdf_binz)

main()
	
