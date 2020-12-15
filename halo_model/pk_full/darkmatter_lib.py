from cosmosis.datablock import names, option_section
import sys
import numpy as np
from scipy.interpolate import interp1d, interp2d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import hankel
from scipy.integrate import quad, simps, trapz
from scipy.special import legendre, sici, binom
import math


# Concentration-mass relations
# At the moment, only Duffy is actually used in the pipeline (first function)
def concentration(mass, z):
        mass_term = 1./((mass/(2.e12))**0.081)
        c_dm = 10.14*mass_term/((1.+z)**1.01)
        return c_dm

def conc_cooraysheth(mass, z, mstar):
        c_dm = (9./(1.+z))*(mass/mstar)**(-0.13)
        return c_dm
        
def conc(m, z, model):
        if model == "Duffy2008":
                mass_term = 1./((m/(2.E+12))**0.081)
                c_dm = 10.14*mass_term/((1.+z)**1.01)
        elif model == "Mandelbaum2008":        
                mass_term = 1./((m/(1.E+14))**0.1)
                c_dm = 5.*mass_term/(1.+z)        
        else:
                print( "Error [conc]: model not specified")        
                #exit(-1)        
        return c_dm


# Virial radius associated with a given mass (r200 in our choice)
def radvir_from_mass(mass, mean_density0):
        return ((3.*mass)/(4.*200.*np.pi*mean_density0))**(1./3.)
        
# Analytic Fourier transform of the Navarro-Frenk-White density profile
def norm_fourier(x, c):
        # Note: x = k*r_scale where r_scale= r200/c
        si, ci = sici((1.+c)*x)
        si2, ci2 = sici(x)
        sinterm = np.sin(x)*( si-si2 ) -np.sin(c*x)/((1+c)*x)
        costerm = np.cos(x)*( ci-ci2 )
        f_c = np.log(1.+c)-c/(1.+c)
        u_fourier = ( sinterm + costerm )/f_c
        return u_fourier

def compute_u_dm(k,mass,z, mean_density0):
        # here the concentration is function of mass for fixed value of z
        c_dm = concentration(mass, z)
        rvir = radvir_from_mass(mass,mean_density0)
        rs = rvir/c_dm
        # compute the analytic fourier-transform of the nfw profile 
        u_dm = norm_fourier(k*rs,c_dm)
        return u_dm        
        
#######        
        
        

        
#######

