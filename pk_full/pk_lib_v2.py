# A new power spectrum module

# NOTE: no truncation nether transition applied!

from cosmosis.datablock import names, option_section
import sys
import numpy as np
from scipy.interpolate import interp1d, interp2d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.integrate import quad, simps, trapz
import math
from itertools import count


from darkmatter_lib import compute_u_dm 


def one_halo_truncation(k_vec):
        return 1.0 #0.5*(1.0+erf((k_vec)-0.01))

def two_halo_truncation(k_vec):
        return 1.0 #0.5*(1.0+(erf(-(k_vec+1.25)))


def unknown_pk(pk_xx):
        raise ValueError("Unknown combination of power spectra. Please, provide a valid one.")
        return

def compute_two_halo_alignment(block, suffix, nz, nk, growth_factor, mean_density0, luminosity_dependence):
        l_term = np.ones(nz)
        l_term_II = np.ones(nz)
        if luminosity_dependence == True:
                l_term = block["ia_lum_scaling", "l_l0_beta"]
                l_term_II = block["ia_lum_scaling", "l_l0_beta2"]
        # linear alignment coefficients
        C1 = 5.e-14
        # load the 2h amplitude of the alignment signal from the value file. Note that gamma_2h is often referred as A_h in literature
        gamma_2h = block["ia_parameters" + suffix, "gamma_2h"]
        alignment_amplitude_2h = np.empty([nz,nk])
        alignment_amplitude_2h_II = np.empty([nz,nk])
        for jz in range (0,nz):
                alignment_amplitude_2h[jz] = -l_term[jz]*(gamma_2h*C1*mean_density0/growth_factor[jz])
                alignment_amplitude_2h_II[jz] = l_term_II[jz]*(gamma_2h*C1*mean_density0/growth_factor[jz])**2
        return alignment_amplitude_2h, alignment_amplitude_2h_II

def compute_central_galaxy_factor(Ncen, numdenscen, f_c):
        return f_c*Ncen/numdenscen

def compute_satellite_galaxy_factor(Nsat, numdenssat, f_s, u_gal):
        return f_s*Nsat*u_gal/numdenssat
        
def compute_matter_factor(mass, mean_density0, u_dm):
        return (mass/mean_density0)*u_dm

def compute_satellite_galaxy_alignment_factor(Nsat, numdenssat, f_s, wkm_sat):
        return f_s*Nsat*wkm_sat/numdenssat

        
def compute_1h_term(factor_1, factor_2, mass, dn_dlnm_z):
        integrand = factor_1*factor_2*dn_dlnm_z/mass
        sum_1h = simps(integrand, mass)
        return sum_1h

def compute_Ig_term(factor_1, mass, dn_dlnm_z, b_m):
        integrand = factor_1*b_m*dn_dlnm_z/mass
        I_g= simps(integrand, mass)
        return I_g
        
def compute_Im_term(mass, u_dm, b_dm, dn_dlnm, mean_density0):
        integrand_m1 = b_dm*dn_dlnm*(1./mean_density0)
        integrand_m2 = b_dm*dn_dlnm*u_dm*(1./mean_density0)
        I_m1 = 1. - simps(integrand_m1, mass)
        I_m2 = simps(integrand_m2, mass)
        I_m = I_m1+I_m2
        return I_m        
        
def compute_2h_term(plin, I1, I2):
        return plin*I1*I2
        
def compute_u_dm_grid(k_vec, mass, z_vec, mean_density0, nz, nk, nmass):
        u_dm = np.empty([nz,nk,nmass])
        for jz in range(0,nz):
                for ik in range(0,nk):
                        u_dm[jz,ik,:] = np.abs(compute_u_dm(k_vec[ik],mass,z_vec[jz], mean_density0))
        return u_dm

# matter
        
def prepare_matter_factor_grid(mass, mean_density0, u_dm, nz, nk, nmass):
        m_factor = np.empty([nz,nk,nmass])        
        for jz in range(0,nz):
                for ik in range(0,nk):        
                        m_factor[jz,ik,:] = compute_matter_factor(mass, mean_density0, u_dm[jz,ik,:])
        return m_factor

def prepare_Im_term(mass, u_dm, b_dm, dn_dlnm, mean_density0, nz, nk):
        I_m_term= np.empty([nz,nk])        
        for jz in range(0,nz):
                for ik in range(0,nk):        
                        I_m_term[jz,ik] = compute_Im_term(mass, u_dm[jz,ik,:], b_dm[jz], dn_dlnm[jz], mean_density0)
        return I_m_term

'''        
def prepare_Im_term(mass, u_dm, b_dm, dn_dlnm, mean_density0, nz, nk):
        I_m_term= np.empty([nz,nk])        
        for jz in range(0,nz):
                for ik in range(0,nk):        
                        I_m_term[jz] = I_m_term[jz].append(compute_Im_term(mass, u_dm[jz,ik,:], b_dm[jz], dn_dlnm[jz], mean_density0), axis=1)
        return I_m_term        
'''        

        
# satellites        

# clustering
def prepare_satellite_factor_grid(mass, Nsat, numdensat, f_sat, u_gal, nz, nk, nmass):
        s_factor = np.empty([nz,nk,nmass])        
        for jz in range(0,nz):
                for ik in range(0,nk):        
                        s_factor[jz,ik,:] = compute_satellite_galaxy_factor(Nsat[jz,:], numdensat[jz], f_sat[jz], u_gal[jz,ik,:])
        return s_factor

def prepare_Is_term(mass, s_factor, b_m, dn_dlnm, nz, nk):
        I_s_term= np.empty([nz,nk])        
        for jz in range(0,nz):
                for ik in range(0,nk):        
                        I_s_term[jz,ik] = compute_Ig_term(s_factor[jz,ik], mass, dn_dlnm[jz], b_m[jz])
        return I_s_term

# alignment       

# This function returns the grid in z, k and mass for the satellite alignment
#
# f_sat/n_sat N_sat gamma_hat(k,M)
# 
# where gamma_hat(k,M) is the Fourier transform of the density weighted shear, i.e. the radial dependent power law times the NFW profile,
# here computed by the module wkm, while gamma_1h is only the luminosity dependence factor.
def prepare_satellite_alignment_factor_grid(mass, Nsat, numdensat, f_sat, wkm, gamma_1h, nz, nk, nmass):
        s_align_factor = np.empty([nz,nk,nmass])        
        for jz in range(0,nz):
                for ik in range(0,nk):        
                        s_align_factor[jz,ik,:] = compute_satellite_galaxy_alignment_factor(Nsat[jz,:], numdensat[jz], f_sat[jz], wkm[jz,:,ik])
                s_align_factor[jz] *= gamma_1h[jz]
        print ('s_align_factor successfully computed!')
        return s_align_factor
        

# centrals

def prepare_central_factor_grid(mass, Ncen, numdencen, f_cen, nz, nmass):
        c_factor = np.array([compute_central_galaxy_factor(Ncen_z, numdencen_z, f_cen_z) for Ncen_z, numdencen_z, f_cen_z in zip(Ncen, numdencen, f_cen)])
        #print 'c_factor.shape = ', c_factor.shape
        return c_factor

def prepare_Ic_term(mass, c_factor, b_m, dn_dlnm, nz, nk):
        I_c_term= np.empty([nz,nk])        
        for jz in range(0,nz):
                I_c_term[jz] = compute_Ig_term(c_factor[jz], mass, dn_dlnm[jz], b_m[jz])
        return I_c_term
        
# ---- POWER SPECTRA ----#

# matter-matter
def compute_p_mm_new(block, k_vec, plin, z_vec, mass, dn_dlnm, m_factor, I_m_term, nz, nk):
        # prepare the 1h term
        pk_mm_1h = np.empty([nz,nk])
        for jz in range(0,nz):
                for ik in range(0,nk):
                        pk_mm_1h[jz,ik] = compute_1h_term(m_factor[jz,ik,:], m_factor[jz,ik,:], mass, dn_dlnm[jz]) 
        # compute the 2h term
        #one_halo_trunc = (1.0+erf(k_vec-0.001))
        #one_halo_trunc = 1.-0.5*(1.+erf(k_vec-2.5)/0.25)
        #two_halo_trunc = 0.5*(1.0+erf(-(k_vec+1.25)))
        #pk_mm_1h = pk_mm_1h #*one_halo_trunc                
        pk_mm_2h = compute_2h_term(plin, I_m_term, I_m_term)#*two_halo_trunc
        pk_mm_tot = pk_mm_1h + pk_mm_2h
        # save in the datablock
        block.put_grid("matter_1h_power", "z", z_vec, "k_h", k_vec, "p_k", pk_mm_1h) 
        block.put_grid("matter_2h_power", "z", z_vec, "k_h", k_vec, "p_k", pk_mm_2h) 
        
        #block.put_double_array_1d("matter_power_nl", "z", z_vec) 
        #block.put_double_array_1d("matter_power_nl", "k_h", k_vec) 
        #block.put_double_array_2d("matter_power_nl", "p_k", pk_mm_tot) 

        
        block.put_grid("matter_power", "z", z_vec, "k_h", k_vec, "p_k", pk_mm_tot) 
        #block.put_grid("matter_power_nl", "k_h", k_vec, "z", z_vec, "p_k", pk_mm_tot) 

        print( 'p_mm succesfully computed')
        return 

# galaxy-galaxy power spectrum
def compute_p_nn(block, k_vec, plin, z_vec, mass, dn_dlnm, c_factor, s_factor, I_c_term, I_s_term, nz, nk):
        #
        # p_tot = p_cs_1h + p_ss_1h + p_cs_2h + p_cc_2h
        #
        # prepare the 1h term
        pk_cs_1h = np.empty([nz,nk])
        pk_ss_1h = np.empty([nz,nk])
        for jz in range(0,nz):
                for ik in range(0,nk):
                        pk_cs_1h[jz,ik] = compute_1h_term(c_factor[jz], s_factor[jz,ik], mass, dn_dlnm[jz]) 
                        pk_ss_1h[jz,ik] = compute_1h_term(s_factor[jz,ik,:], s_factor[jz,ik,:], mass, dn_dlnm[jz])
                mask_large_scales = k_vec < 0.01
                pk_cs_1h[jz][mask_large_scales] = 0.0
                pk_ss_1h[jz][mask_large_scales] = 0.0
        pk_cs_2h = compute_2h_term(plin, I_c_term, I_s_term)#*two_halo_trunc
        pk_cc_2h = compute_2h_term(plin, I_c_term, I_c_term)#*two_halo_trunc
        pk_ss_2h = compute_2h_term(plin, I_s_term, I_s_term)#*two_halo_trunc
        pk_tot = 2. * pk_cs_1h + pk_ss_1h + pk_cc_2h + pk_ss_2h + 2.*pk_cs_2h 
        # save in the datablock
        #block.put_grid("galaxy_cc_power_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_cc_2h) 
        #block.put_grid("galaxy_ss_power_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_ss_2h) 
        #block.put_grid("galaxy_cs_power_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_cs_2h) 
        block.put_grid("galaxy_cs_power_1h", "z", z_vec, "k_h", k_vec, "p_k", pk_cs_1h) 
        block.put_grid("galaxy_ss_power_1h", "z", z_vec, "k_h", k_vec, "p_k", pk_ss_1h) 
        #block.put_grid("galaxy_power", "z", z_vec, "k_h", k_vec, "p_k", pk_tot) 
        #block.put_grid("galaxy_power", "k_h", k_vec, "z", z_vec, "p_k", pk_tot) 

        # galaxy linear bias
        #block.put_grid("galaxy_linear_bias_cen", "z", z_vec, "k_h", k_vec, "galaxybiascentralsred", I_c_term)
        #block.put_grid("galaxy_linear_bias_sat", "z", z_vec, "k_h", k_vec, "galaxybiassatellitesred", I_s_term)
        #block.put_grid("galaxy_linear_bias_cenxsat", "z", z_vec, "k_h", k_vec, "galaxybiascentralsxsatellitesred", np.sqrt(I_s_term*I_c_term))
        #block.put_grid("galaxy_linear_bias_tot", "z", z_vec, "k_h", k_vec, "galaxybiastotalred", np.sqrt(I_c_term**2.+I_s_term**2.+2.*I_s_term*I_c_term))
        galaxy_linear_bias = np.sqrt(I_c_term**2.+I_s_term**2.+2.*I_s_term*I_c_term)
        
        print( 'p_nn succesfully computed')
        return pk_cs_1h+pk_ss_1h, pk_cc_2h+pk_ss_2h+pk_cs_2h, pk_tot, galaxy_linear_bias

# galaxy-matter power spectrum
def compute_p_xgG(block, k_vec, plin, z_vec, mass, dn_dlnm, c_factor, s_factor, m_factor, I_c_term, I_s_term, I_m_term, nz, nk):
        #
        # p_tot = p_cm_1h + p_sm_1h + p_cm_2h + p_cm_2h
        #
        # prepare the 1h term
        pk_cm_1h = np.empty([nz,nk])
        pk_sm_1h = np.empty([nz,nk])
        for jz in range(0,nz):
                for ik in range(0,nk):
                        pk_cm_1h[jz,ik] = compute_1h_term(c_factor[jz], m_factor[jz,ik], mass, dn_dlnm[jz]) 
                        pk_sm_1h[jz,ik] = compute_1h_term(s_factor[jz,ik], m_factor[jz,ik], mass, dn_dlnm[jz]) 
        #one_halo_trunc = 0.5*(1.0+erf(k_vec-0.001))
        #two_halo_trunc = 0.5*(1.0+erf(-(k_vec+1.25)))
        #pk_cs_1h = pk_cs_1h#*one_halo_trunc
        #pk_ss_1h = pk_ss_1h#*one_halo_trunc
        pk_cm_2h = compute_2h_term(plin, I_c_term, I_m_term)#*two_halo_trunc
        pk_sm_2h = compute_2h_term(plin, I_s_term, I_m_term)#*two_halo_trunc
        pk_tot = pk_cm_1h + pk_sm_1h + pk_cm_2h + pk_sm_2h
        block.put_grid("matter_cm_galaxy_power_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_cm_2h) 
        block.put_grid("matter_sm_galaxy_power_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_sm_2h) 
        block.put_grid("matter_cm_galaxy_power_1h", "z", z_vec, "k_h", k_vec, "p_k", pk_cm_1h) 
        block.put_grid("matter_sm_galaxy_power_1h", "z", z_vec, "k_h", k_vec, "p_k", pk_sm_1h) 
        block.put_grid("matter_galaxy_power", "z", z_vec, "k_h", k_vec, "p_k", pk_tot) 
        print( 'p_xgG succesfully computed')
        return 
        
        


#################################################
#                                                                                                #
#                INTRINSIC ALIGNMENT POWER SPECTRA                #
#                                                                                                #        
#################################################


# galaxy-matter power spectrum
def compute_p_xGI(block, k_vec, p_eff, z_vec, mass, dn_dlnm, m_factor, s_align_factor, alignment_amplitude_2h, nz, nk, f_gal):
        #
        # p_tot = p_sm_GI_1h + f_cen*p_cm_GI_2h + O(any other combination)
        #
        # prepare the 1h term
        pk_sm_1h = np.empty([nz,nk])
        for jz in range(0,nz):
                for ik in range(0,nk):
                        # here we are assuming the convention where the radial alignment is negative
                        pk_sm_1h[jz,ik] = - compute_1h_term(m_factor[jz,ik,:], s_align_factor[jz,ik,:], mass, dn_dlnm[jz]) 
        # this is simply the Linear Alignment Model
        pk_cm_2h = compute_p_xGI_two_halo(block, k_vec, p_eff, z_vec, nz, f_gal, alignment_amplitude_2h)
        pk_tot = pk_sm_1h*one_halo_truncation(k_vec) + pk_cm_2h*two_halo_truncation(k_vec)
        # save in the datablock
        #block.put_grid("matter_intrinsic_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_cm_2h) 
        #block.put_grid("matter_intrinsic_1h", "z", z_vec, "k_h", k_vec, "p_k", pk_sm_1h) 
        #block.put_grid("matter_intrinsic_power", "z", z_vec, "k_h", k_vec, "p_k", pk_tot) 
        print( 'p_xGI succesfully computed')
        return pk_sm_1h, pk_cm_2h, pk_tot



# galaxy-intrinsic power spectrum
def compute_p_gI(block, k_vec, p_eff, z_vec, mass, dn_dlnm, c_factor, s_align_factor, I_c_term, alignment_amplitude_2h, nz, nk):
        #
        # p_tot = p_cs_gI_1h + (2?)*p_cc_gI_2h + O(p_ss_gI_1h) + O(p_cs_gI_2h)
        #
        # prepare the 1h term
        pk_cs_1h = np.empty([nz,nk])
        for jz in range(0,nz):
                for ik in range(0,nk):
                        pk_cs_1h[jz,ik] = compute_1h_term(c_factor[jz], s_align_factor[jz,ik,:], mass, dn_dlnm[jz]) 
        #print "did it work?"
        # this is simply the Linear Alignment Model
        # NOTE: that here we are assuming the linear bias to be caused by the central galaxies only
        # if we want to use the bias of the entire sample, we can obtain it as:
        # bias = np.sqrt(Ic**2 + Is**2 + (2.*Ic*Is) [see the notes on the spiral notebook] 
        align_2h = np.empty(I_c_term.shape)
        for jz in range (0,nz):
                align_2h[jz] = alignment_amplitude_2h[jz]
        #b_g = np.sqrt(I_c_term**2.+I_s_term**2.+2.*I_s_term*I_c_term)
        #pk_tot_2h =  b_g*(p_eff*align_2h)
        pk_cc_2h =  compute_2h_term(p_eff, I_c_term, align_2h)
        pk_tot = pk_cs_1h + pk_cc_2h
        # save in the datablock
        block.put_grid("galaxy_cc_intrinsic_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_cc_2h) 
        block.put_grid("galaxy_cs_intrinsic_1h", "z", z_vec, "k_h", k_vec, "p_k", pk_cs_1h) 
        block.put_grid("galaxy_intrinsic_power", "z", z_vec, "k_h", k_vec, "p_k", pk_tot) 
        print( 'p_gI succesfully computed')
        return 
        
# intrinsic-intrinsic power spectrum
def compute_p_II(block, k_vec, p_eff, z_vec, mass, dn_dlnm, s_align_factor, alignment_amplitude_2h_II, nz, nk, f_gal):
        #
        # p_tot = p_ss_II_1h + p_cc_II_2h + O(p_sc_II_1h) + O(p_cs_II_2h)
        #
        # prepare the 1h term
        pk_ss_1h = np.empty([nz,nk])
        for jz in range(0,nz):
                for ik in range(0,nk):
                        pk_ss_1h[jz,ik] = compute_1h_term(s_align_factor[jz,ik,:], s_align_factor[jz,ik,:], mass, dn_dlnm[jz]) # - compute_1h_term(s_align_factor[jz,ik,:], 1.0, mass, dn_dlnm[jz]) )
        # This is simply the Linear Alignment Model weighted by the central galaxy fraction
        pk_tot_2h = compute_p_II_two_halo(block, k_vec, p_eff, z_vec, nz, f_gal, alignment_amplitude_2h_II)
        pk_tot = np.empty([nz,nk])
        for jz in range(0,nz):
                pk_tot[jz] = pk_ss_1h[jz] + pk_tot_2h[jz]
        # save in the datablock
        #block.put_grid("galaxy_cc_intrinsic_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_cc_2h) 
        #block.put_grid("intrinsic_power_2h", "z", z_vec, "k_h", k_vec, "p_k", pk_tot_2h) 
        #block.put_grid("intrinsic_power_1h", "z", z_vec, "k_h", k_vec, "p_k", pk_ss_1h) 
        #block.put_grid("intrinsic_power", "z", z_vec, "k_h", k_vec, "p_k", pk_tot) 
        print( 'p_II succesfully computed')
        return pk_ss_1h, pk_tot_2h, pk_tot
                
        
############### TWO HALO ONLY ###################

# galaxy-galaxy power spectrum
def compute_p_nn_two_halo(block, k_vec, plin, z_vec, bg):
        #
        # p_tot = b_g**2 * p_lin
        #
        pk_tot = np.zeros([len(z_vec),len(k_vec)])
        for jz in range(len(z_vec)):
                pk_tot[jz] = bg[jz]**2.*plin[jz] 
        return pk_tot

# galaxy-matter power spectrum
def compute_p_xgG_two_halo(block, k_vec, plin, z_vec, bg):
        #
        # p_tot = bg * plin
        #
        pk_tot = np.zeros([len(z_vec),len(k_vec)])
        for jz in range(len(z_vec)):
                pk_tot[jz] = bg[jz]*plin[jz]         
        return pk_tot

# galaxy-matter power spectrum
def compute_p_xGI_two_halo(block, k_vec, p_eff, z_vec, nz, f_gal, alignment_amplitude_2h):
        #
        # p_tot = p_NLA
        #
        # this is simply the Linear (or Nonlinear) Alignment Model, weighted by the central galaxy fraction
        pk_tot = np.zeros([len(z_vec),len(k_vec)])
        for jz in range (0,nz):
                pk_tot[jz] = f_gal[jz]*alignment_amplitude_2h[jz]*p_eff[jz]
        return pk_tot


# galaxy-intrinsic power spectrum
def compute_p_gI_two_halo(block, k_vec, p_eff, z_vec, nz, f_gal, alignment_amplitude_2h, bg):
        #
        # p_tot = bg * p_NLA
        #
        pk_tot = np.zeros([len(z_vec),len(k_vec)])
        for jz in range (0,nz):
                pk_tot[jz] = f_gal[jz]*bg[jz]*alignment_amplitude_2h[jz]*p_eff[jz]        
        return pk_tot
        
# galaxy-intrinsic power spectrum
def compute_p_II_two_halo(block, k_vec, p_eff, z_vec, nz, f_gal, alignment_amplitude_2h_II):
        pk_tot = np.zeros([len(z_vec),len(k_vec)])
        for jz in range (0,nz):
                pk_tot[jz] = (f_gal[jz]**2.)*p_eff[jz]*alignment_amplitude_2h_II[jz]
        return pk_tot
        
        
        
#################################################
