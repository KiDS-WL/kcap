# Compute the HMF through the Steven Murray python package hmf
# At the moment this module is ONLY designed to work in this particular pipeline
# To wrap it properly in CosmoSIS more work is needed. It also requires a
# permission from the authors -> I haven't contact them yet, so this must be
# considered only for a private usage.

# The module also includes an option to compute the halo bias from Tinker+10.



from cosmosis.datablock import names, option_section
import sys
import numpy as np
from scipy.interpolate import interp1d, interp2d
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import simps
import hmf
from hmf import MassFunction


def tinker_bias(nu, Delta=200., delta_c=1.686):
	# Table 2, Tinker+2010
	y = np.log10(Delta)
	expvar = np.exp(-(4./y)**4.)
	A = 1.+0.24*y*expvar
	a = 0.44*y-0.88
	B = 0.183
	b = 1.5
	C = 0.019+0.107*y+0.19*expvar
	c = 2.4
	# equation 6
	bias = 1.-A*(nu**a)/(nu**a+delta_c**a) + B*nu**b + C*nu**c
	return bias

def normalise_hbf_nu(b_nu, f_nu, nu):
	norm = simps(b_nu*f_nu,nu)
	return b_nu/norm

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.
	
	log_mass_min = options[option_section, "log_mass_min"]
	log_mass_max = options[option_section, "log_mass_max"]
	nmass = options[option_section, "nmass"]
	
	zmin = options[option_section, "zmin"]
	zmax = options[option_section, "zmax"]
	nz = options[option_section, "nz"]
	z_vec = np.linspace(zmin, zmax, nz)
		
	#dlog10m = (log_mass_max-log_mass_min)/nmass
	
	#hmf_model = options[option_section, "hmf_model"]
	halo_bias_option = options[option_section, "do_halo_bias"]

	return log_mass_min, log_mass_max, nmass, z_vec, nz, halo_bias_option #dlog10m,



def execute(block, config):
    #This function is called every time you have a new sample of cosmological and other parameters.
    #It is the main workhorse of the code. The block contains the parameters and results of any 
    #earlier modules, and the config is what we loaded earlier.
	
	log_mass_min, log_mass_max, nmass, z_vec, nz, halo_bias_option = config #dlog10m, 

	# Load cosmological parameters
	cosmology = FlatLambdaCDM(H0=block[cosmo, "hubble"], Om0=block[cosmo, "omega_m"], Ob0=block[cosmo, "omega_b"])
	ns = block[cosmo, "n_s"]
	
	#--------------------------------------#
	# read sigma_8 for the given cosmology
	#--------------------------------------#

	# Note that CAMB does not return the sigma_8 at z=0, as it might seem from the documentation, but sigma_8(z), 
	# so it can't be read from the datablock.
	# This is a problem when fitting for cosmological parameters.
	# A possible solution is to use the rescale_sigma8_module.
	# Be careful that there have been discussions around this module, so I do not know how trustable it is.

	sigma_8 = block[cosmo, "sigma_8"]
	print 'sigma_8 = ', sigma_8
	
	nmass_hmf = 500

	dndlnmh = np.empty([nz, nmass_hmf])
	nu = np.empty([nz,nmass_hmf])
	# Note that fsigma, as function of sigma, is a one dim array -> the dependence in z is absorbed in sigma
	fsigma = np.empty([nz, nmass_hmf])
	fnu = np.empty([nz, nmass_hmf])
	
	# use a large range in masses to properly normalise the halo bias function
	dlog10m = (16.-2.)/500.

	for jz in range(0,nz):
		mf = MassFunction(z= z_vec[jz], Mmin=2., Mmax=16., dlog10m=dlog10m, \
				cosmo_model=cosmology, sigma_8 = sigma_8, n = ns, hmf_model='Tinker10')
		nu[jz] = mf.nu**0.5
		dndlnmh[jz] = mf.dndlnm
		fsigma[jz] = mf.fsigma
		fnu[jz] = fsigma[jz]/nu[jz]
	
	mass = np.logspace(log_mass_min, log_mass_max, nmass)
	f_interp_hmf = interp2d(mf.m, z_vec, dndlnmh)
	hmf_interp = f_interp_hmf(mass, z_vec)
	block.put_grid("hmf", "z", z_vec, "m_h", mass, "dndlnmh", hmf_interp)
	

	#--------------------------------------#	
	# HALO BIAS
	#--------------------------------------#

	Delta = 200.
	delta_c = 1.686 
	
	if halo_bias_option:
		#b_nu = np.empty([nz,len(mf.m)])
		hb_normalised = np.empty([nz,len(mf.m)])

		for jz in range(0,nz):
			b_nu = tinker_bias(nu[jz], Delta, delta_c)
			hb_normalised[jz] = normalise_hbf_nu(b_nu, fnu[jz], nu[jz])
		f_interp_hb = interp2d(mf.m, z_vec, hb_normalised)
		hb_interp = f_interp_hb(mass, z_vec)
		block.put_grid("halobias", "z", z_vec, "m_h", mass, "b_hb", hb_interp)


	#if save_all:
	#	block.put_grid("hmf", "z", z_vec, "m_h", mf.m, "nu", nu)
	#	block.put_grid("hmf", "z", z_vec, "m_h", mf.m, "fsigma", fsigma)
	#	block.put_grid("hmf", "z", z_vec, "m_h", mf.m, "fnu", fnu)

	return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
