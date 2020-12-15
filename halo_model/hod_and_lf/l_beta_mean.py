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


from cosmosis.datablock import names, option_section
import sys
import re
from os import listdir
import numpy as np
import clf_lib as lf
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import simps

import time

cosmo = names.cosmological_parameters


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
	
def mean_L_L0_to_beta(loglum, l0, beta):
    lum = 10.**loglum
    beta2 = 2*beta
    pdf, lum_bins = np.histogram(lum, bins=800)
    pdf_test, lum_bins = np.histogram(lum, bins=800, density=True)
    #print (np.sum(pdf_test*np.diff(lum_bins)))
    dbin = (lum_bins[-1]-lum_bins[0])/800.
    bincen = lum_bins[0:-1]+0.5*dbin
    norm_pdf = simps(pdf, bincen)
    l_beta_mean = simps(pdf*(bincen/l0)**beta, bincen)/norm_pdf
    #print ('l_beta_mean = ', l_beta_mean)
    #print ('l_beta_mean_pdf = ', np.sum(pdf_test*(bincen/l0)**beta)*np.diff(lum_bins))
    l_beta2_mean = simps(pdf*(bincen/l0)**beta2, bincen)/norm_pdf
    return l_beta_mean, l_beta2_mean

def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.
	
	# general version, to be implemented after I finish with this project
	# split the data into the zbins

	directory_name = options[option_section, "path_to_lum_files"] 
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
		
	return zbin, loglum



def execute(block, config):

	#loglum_bin1, loglum_bin2, loglum_bin3, loglum_bin4, loglum_bin5, loglum_bin6, nz = config
	
	zbin, loglum = config
	
	nz = len(zbin)
	
	l0 = block["ia_luminosity_scaling", "l_0"]
	beta_l = block["ia_luminosity_scaling", "beta_l"]
	
	mean_lscaling = np.empty(nz)
	mean_lscaling_beta2 = np.empty(nz)
	

	for i in range(0,nz):
		z_i = 'z%s' %str(i+1)
		mean_lscaling[i], mean_lscaling_beta2[i] = mean_L_L0_to_beta(loglum[z_i], l0, beta_l)
	

	block.put_double_array_1d("ia_lum_scaling", "l_l0_beta", mean_lscaling)
	block.put_double_array_1d("ia_lum_scaling", "l_l0_beta2", mean_lscaling_beta2)

	return 0
	

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
	
	
