from cosmosis.datablock import names, option_section
import sys
import numpy as np
import clf_lib as lf
from scipy.interpolate import interp1d, interp2d, RegularGridInterpolator
from scipy.integrate import simps

import time


# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def setup(options):
	nz_tomographic = options[option_section, "nz_tomographic"] 
	nz_output_section = options[option_section, "nz_output_section"] 
	spec_option = options[option_section, "spec_option"]
	print (spec_option)
	z_median = -99.
	if spec_option == True:
		z_median = options[option_section, "z_median"]
	return nz_tomographic, nz_output_section, spec_option, z_median
	

def execute(block, config):
	# read the inputs from the block (if any)
	nz_tomographic, nz_output_section, spec_option, z_median = config
	# read lf from the block:
	z_lf = block["luminosity_function", "z"]
	lum = block["luminosity_function", "lum"]
	lf_l = block["luminosity_function", "lf_l"]
	# intepolate
	interp_lf = interp2d(lum, z_lf, lf_l) 
	if spec_option == True:
		lf_l_photz = interp_lf(lum, z_median)
		block.put_double_array_1d("photometric_luminosity_function", "lum", lum)
		block.put_double_array_1d("photometric_luminosity_function", "lf_l_photz", lf_l_photz)

	else:
		# read nz from the block:
		bin_z = {}
		bin = [] # this is for the output format
		for i in range(1,nz_tomographic+1):
			bin_z["{0}".format(i-1)] = block[nz_output_section, "bin_%d" %i]
			#if np.isnan(bin_z['%d' %(i-1)]):
			#	print ("Error! Invalid pdf!")
			bin = np.append(bin, i)
		z_array = block[nz_output_section, "z"]
		lf_sampled = interp_lf(lum, z_array)
	
		l_phot_z = np.empty([nz_tomographic, len(lum)])
		# integrate
		for jz in range(0,nz_tomographic):
			for il in range(0,len(lum)):
				l_phot_z[jz, il] = simps(lf_sampled[:,il]*bin_z['%d' %jz], z_array)
		# save the results in the datablock
		block.put_grid("photometric_luminosity_function", "bin", bin, "lum", lum, "lf_l_photz", l_phot_z  )
	return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
