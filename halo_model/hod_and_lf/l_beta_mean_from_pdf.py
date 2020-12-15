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
	data = np.loadtxt(file_name, dtype=np.float).T
	return data


def mean_L_L0_to_beta(xlum, pdf, l0, beta):
    l_beta_mean = simps(pdf*(xlum/l0)**beta, xlum)
    l_beta2_mean = simps(pdf*(xlum/l0)**(2.*beta), xlum)
    return l_beta_mean, l_beta2_mean

def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.
	
	# general version, to be implemented after I finish with this project
	# split the data into the zbins

	filename = options[option_section, "path_to_lum_file"] 
	galaxy_type_option = options[option_section, "galaxy_type"]
	
	if galaxy_type_option == "centrals":
		galaxy_type = 0
	elif galaxy_type_option == "satellites":
		galaxy_type = 1
	else:
		raise ValueError("Please, select galaxy type. Options are: 'centrals' or 'satellites'.")
	
	data = load_data(filename)
	
	lum = data[0]
	nlum = len(lum)
	nz = np.size(data, axis=0) -1
	print ('nz = %d' %nz)
	lum_pdf_z = data[1:,:]
	
		
	name = options.get_string(option_section, "name", default="").lower()
	if name:
		suffix = "_" + name
	else:
		suffix = ""
	
	return lum, lum_pdf_z, nz, nlum, galaxy_type, suffix



def execute(block, config):

	lum, lum_pdf_z, nz, nlum, galaxy_type, suffix = config
		
	if galaxy_type==0:		
		l0 = block["ia_luminosity_scaling", "l_0"]
		beta_l = block["ia_luminosity_scaling", "beta_l"]	
		mean_lscaling = np.empty(nz)
		mean_lscaling_beta2 = np.empty(nz)
		for i in range(0,nz):
			mean_lscaling[i], mean_lscaling_beta2[i] = mean_L_L0_to_beta(lum, lum_pdf_z[i], l0, beta_l)
		block.put_double_array_1d("ia_lum_scaling", "l_l0_beta", mean_lscaling)
		block.put_double_array_1d("ia_lum_scaling", "l_l0_beta2", mean_lscaling_beta2)
		
	if galaxy_type==1:		
		l0 = block["ia_satellite_luminosity_scaling" + suffix, "l_0"]
		beta_l = block["ia_satellite_luminosity_scaling" + suffix, "beta_l"]	
		mean_lscaling = np.empty(nz)
		# this will be introduced in case an II term will be implemented
		# mean_lscaling_beta2 = np.empty(nz)
		for i in range(0,nz):
			mean_lscaling[i], _ = mean_L_L0_to_beta(lum, lum_pdf_z[i], l0, beta_l)
		block.put_double_array_1d("ia_satellite_lum_scaling" + suffix, "l_l0_beta", mean_lscaling)
		#block.put_double_array_1d("ia_lum_scaling", "l_l0_beta2", mean_lscaling_beta2)

	return 0
	

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
	
	
