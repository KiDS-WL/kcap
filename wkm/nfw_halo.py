from cosmosis.datablock import names, option_section
import sys
import numpy as np

import time

import math


from uell_radial_dependent_alignment_lib import minimal_cosmo, radvir

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.

cosmo = names.cosmological_parameters


#--------------------------------------------------------------------------------#	

def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.
	
	log_mass_min = options[option_section, "log_mass_min"]
	log_mass_max = options[option_section, "log_mass_max"]
	nmass = options[option_section, "nmass"]
	#log-spaced mass in units of M_sun/h
	mass = np.logspace(log_mass_min, log_mass_max, nmass)
	
	zmin = options[option_section, "zmin"]
	zmax = options[option_section, "zmax"]
	nz = options[option_section, "nz"]
	z_vec = np.linspace(zmin, zmax, nz)
	print(z_vec)
	 	 	 
	return z_vec, nz, mass, nmass
	
def execute(block, config):
    #This function is called every time you have a new sample of cosmological and other parameters.
    #It is the main workhorse of the code. The block contains the parameters and results of any 
    #earlier modules, and the config is what we loaded earlier.
    
	z, nz, mass, nmass = config
	start_time = time.time()


	#---- Comology ----#
	# Load cosmological parameters
	omega_m = block[cosmo, "omega_m"]

	#concentration -> to be replaced by a module
	c = np.empty([nz,nmass])
	mass_term = 10.14/((mass/(2.e+12))**0.081)
	for j in range (0,nz):
		for i in range (0,nmass):
			c[j,i]= mass_term[i]/((1.+z[j])**1.01)

	# compute the virial radius and the scale radius associated with a halo of mass M
	rho_m = minimal_cosmo(omega_m)
	rho_halo = 200.*rho_m
	rvir = radvir(mass, rho_halo)
	r_s = np.empty([nz, nmass])
	for j in range(0,nz):
		r_s[j] = rvir/c[j]
		
	#print c.shape
	#print rvir.shape
	#print r_s.shape
	#print z
	#print mass
		
	block.put_grid("concentration", "z", z, "m_h", mass, "c", c)
	#print "ok"
	block.put_grid("nfw_scale_radius", "z", z, "m_h", mass, "rs", r_s)
	#print "ok"
	block.put_double_array_1d("virial_radius", "rvir", rvir)
	#print "ok"

	return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
