from cosmosis.datablock import names, option_section
import sys
import numpy as np
#from os import listdir
#import re

from cosmosis.datablock import names, option_section
#import lf_lib_simps as lf
#from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator
#from itertools import count, izip

from scipy import interp
import matplotlib.pyplot as plt

import time


# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def	extrapolate_z(z_ext, z_vec, pk, nk):
	nz_ext = len(z_ext)
	pk_extz = np.empty([nz_ext, nk])
	for ik in range(0,nk):
		pk_kfixed = pk[:,ik]
		pk_extz[:,ik] = interp(z_ext, z_vec, pk_kfixed)
	return pk_extz

def	extrapolate_k(k_ext, k, pk, nz):
	nk_ext = len(k_ext)
	pk_extk = np.empty([nz, nk_ext])
	for jz in range(0,nz):
		pk_zfixed = pk[jz,:]
		pk_extk[jz,:] = interp(k_ext, k, pk_zfixed)
	return pk_extk

def add_red_and_blue_power(block, f_red, power_section, zmax):
		k = block[power_section+"_red", "k_h"]
		z = block[power_section+"_red", "z"]
		nz = len(z)
		nk = len(k)
		pk_tot = np.zeros([nz,nk])
		pk_red = block[power_section+"_red", "p_k"]
		pk_blue = block[power_section+"_blue", "p_k"]
		for jz in range(nz):
			pk_tot[jz] = f_red[jz]*pk_red[jz] + (1.-f_red[jz])*pk_blue[jz]
			
		# extrapolate
		nz_ext = 100
		z_ext = np.linspace(0., zmax, nz_ext)
		nk_ext = 2580
		k_ext = np.logspace(-4, 4, nk_ext)
		pk_tot_ext_z = extrapolate_z(z_ext, z, pk_tot, nk)
		pk_tot_ext = extrapolate_k(k_ext, k, pk_tot_ext_z, nz_ext)
		#for i in range(0,nz_ext):
		#	plt.loglog(k_ext, np.abs(pk_tot_ext[i]))
		#plt.show()
		block.put_grid(power_section, "z", z_ext, "k_h", k_ext, "p_k", pk_tot_ext)
		
#--------------------------------------------------------------------------------#	

def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.
		
	f_red_file = options[option_section, "f_red_file"]
	f_red = np.loadtxt(f_red_file, unpack=True, usecols=(1,))
	print f_red
	# clustering
	p_nn_option = options[option_section, "do_p_nn"]
	# galaxy lensing
	p_xgG_option = options[option_section, "do_p_xgG"]
	# intrinsic alignment
	p_xGI_option = options[option_section, "do_p_xGI"]
	p_II_option = options[option_section, "do_p_II"]
	p_gI_option = options[option_section, "do_p_gI"]

	zmax =  options[option_section, "zmax"]
			
	return f_red, p_nn_option, p_xgG_option, p_xGI_option, p_II_option, p_gI_option, zmax
	

def execute(block, config):
    #This function is called every time you have a new sample of cosmological and other parameters.
    #It is the main workhorse of the code. The block contains the parameters and results of any 
    #earlier modules, and the config is what we loaded earlier.
	
	f_red, p_nn_option, p_xgG_option, p_xGI_option, p_II_option, p_gI_option, zmax = config

	# TODO: add an error message if len(f_red) != len(z) of the power spectra
	#import pdb; pdb.set_trace()	
	print f_red.shape
	if p_nn_option == True:
		add_red_and_blue_power(block, f_red, "galaxy_power", zmax)	
	if p_xgG_option == True:
		add_red_and_blue_power(block, f_red, "matter_galaxy_power", zmax)	
	if p_xGI_option == True:
		add_red_and_blue_power(block, f_red, "matter_intrinsic_power", zmax)	
	if p_xGI_option == True:
		add_red_and_blue_power(block, f_red, "matter_intrinsic_power", zmax)
	if p_II_option == True:
		add_red_and_blue_power(block, f_red, "intrinsic_power", zmax)
	if p_gI_option == True:
		add_red_and_blue_power(block, f_red, "galaxy_intrinsic_power", zmax)
				
	return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass


