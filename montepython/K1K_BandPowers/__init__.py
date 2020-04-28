##############################################################
# Likelihood for the KiDS-1000 cosmic shear band powers      #
##############################################################
#
# To be used with the KCAP likelihood package for the fiducial
# KiDS1000 COSEBIs analysis from XXXX et al. 2020
# (arXiv:20XX.YYYY).
# Modified to work with Monte Python by Fabian Koehlinger.
#
# Data available from:
#
# http://kids.strw.leidenuniv.nl/sciencedata.php
#
# ATTENTION:
# This likelihood only produces valid results for \Omega_k = 0,
# i.e. flat cosmologies!
##############################################################

from __future__ import print_function
from montepython.likelihood_class import Likelihood
#import io_mp

import os
import sys
import numpy as np
from scipy import interpolate as itp
from scipy.linalg import cholesky, solve_triangular
from astropy.io import fits
# check explicitly for CosmoSIS install:
try:
    import cosmosis.runtime.module
    import cosmosis.datablock
except:
    raise ImportError("Could not import required CosmoSIS modules! \n You can get it from here: " +
                      "\n git clone https://bitbucket.org/tilmantroester/cosmosis.git@kcap#egg=cosmosis-standalone" +
                      "\n \n Follow their README instructions to set it up on your system!")
#from timeit import default_timer as timer

# Python 2.x - 3.x compatibility: Always use more efficient range function
try:
    xrange
except NameError:
    xrange = range

def block_print():
    sys.stdout = open(os.devnull, 'w')

def enable_print():
    sys.stdout = sys.__stdout__

def dict_to_datablock(d={}):
    """
    Convenience function for passing a Python dictionary to a CosmoSIS.datablock
    structure.

    Copy & pasted from kcap/kcap/cosmosis_utils.py

    By Tilman Troester
    """

    b = cosmosis.datablock.DataBlock()
    for section in d.keys():
        for name, value in d[section].items():
            b[section, name] = value

    return b

class K1K_BandPowers(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # Force the cosmological module to store Pk for redshifts up to
        # max(self.z) and for k up to k_max
        self.need_cosmo_arguments(data, {'output': 'mPk'})
        self.need_cosmo_arguments(data, {'P_k_max_h/Mpc': self.k_max_h_by_Mpc})

        ## Compute non-linear power spectrum if requested
        # it seems like HMcode needs the full argument to work...
        if self.method_non_linear_Pk in ['halofit', 'HALOFIT', 'Halofit', 'hmcode', 'Hmcode', 'HMcode', 'HMCODE']:
            self.need_cosmo_arguments(data, {'non linear': self.method_non_linear_Pk})
            print('Using {:} to obtain the non-linear P(k, z)! \n'.format(self.method_non_linear_Pk))
        else:
            print('Only using the linear P(k, z) for ALL calculations \n (check keywords for "method_non_linear_Pk"). \n')

        # TODO: move min_kmax_hmc to data-file?!
        # might not be really necessary; I didn't see a difference in the P(k, z) ratios between
        # HMcode complaining about k_max being too low and not complaining at all...
        if self.method_non_linear_Pk in ['hmcode', 'Hmcode', 'HMcode', 'HMCODE']:
            #self.need_cosmo_arguments(data, {'hmcode_min_k_max': 1000.})
            min_kmax_hmc = 170.
            if self.k_max_h_by_Mpc < min_kmax_hmc:
                self.need_cosmo_arguments(data, {'P_k_max_h/Mpc': min_kmax_hmc})
                #print "Your choice of k_max_h_by_Mpc is too small for HMcode. \n Requested P_k_max_h/Mpc now up to k = {:.2f} h/Mpc \n This does NOT influence the scale above".format(min_kmax_hmc)

        # set up array of ells for Cl integrations:
        self.ells = np.logspace(np.log10(self.ell_min), np.log10(self.ell_max), self.nells)

        # set up KCAP's scale cuts module here:
        # Initialize the scale cuts module from CosmoSIS:
        self.config_scale_cuts = {'scale_cuts': {'data_and_covariance_fits_filename': os.path.join(self.data_directory, self.data_file),
                                                 #'scale_cuts_option': self.scale_cuts_option1,
                                                 'use_stats': 'PeeE',
                                                 'output_section_name': 'likelihood_bp',
                                                 'bandpower_e_cosmic_shear_extension_name': 'PeeE',
                                                 'bandpower_e_cosmic_shear_section_name': 'bandpower'}}

        # for now we only look for these two keywords:
        if hasattr(self, 'keep_ang_PeeE'):
            self.config_scale_cuts['scale_cuts'].update({'keep_ang_PeeE': self.keep_ang_PeeE})
        if hasattr(self, 'cut_pair_PeeE'):
            self.config_scale_cuts['scale_cuts'].update({'cut_pair_PeeE': self.cut_pair_PeeE})

        # import scale_cuts as CosmoSIS module
        self.scale_cuts_module = cosmosis.runtime.module.Module(module_name='scale_cuts',
                                            file_path=os.path.join(self.kcap_directory,
                                            'modules/scale_cuts/scale_cuts.py'))

        # during set up the module stores the cut data vec and covmat in its data
        # attribute
        self.scale_cuts_module.setup(dict_to_datablock(self.config_scale_cuts))

        # this works now:
        self.data_vec = self.scale_cuts_module.data['data']
        #print(self.data_vec.shape)
        covmat = self.scale_cuts_module.data['covariance']
        #print(covmat.shape)

        # Read measurements E_n COSEBIs:
        #self.cosebis_obs = self.__load_legacy_vec_obs()

        # Read covariance matrix
        #covmat = self.__load_legacy_covmat()

        # precompute Cholesky transform for chi^2 calculation:
        self.cholesky_transform = cholesky(covmat, lower=True)

        # Read dn_dz from data FITS file:
        # z_samples, hist_samples = self.__load_legacy_nofz()
        # in the process we also set self.nzbins!
        self.z_samples, self.hist_samples = self.load_data_file()

        # Check if any of the n(z) needs to be shifted in loglkl by D_z{1...n}:
        self.shift_n_z_by_D_z = np.zeros(self.nzbins, 'bool')
        for zbin in xrange(self.nzbins):
            param_name = 'D_z{:}'.format(zbin + 1)
            if param_name in data.mcmc_parameters:
                self.shift_n_z_by_D_z[zbin] = True

        if self.shift_n_z_by_D_z.any():
            # load the correlation matrix of the D_z shifts:
            try:
                fname = os.path.join(self.data_directory, self.filename_corrmat_D_z)
                corrmat_D_z = np.loadtxt(fname)
                print('Loaded correlation matrix for D_z shifts from: \n {:} \n'.format(fname))
                self.L_matrix_D_z = np.linalg.cholesky(corrmat_D_z)
            except:
                print('Could not load correlation matrix of D_z shifts, hence treating them as independent! \n')
                self.L_matrix_D_z = np.eye(self.nzbins)

        # prevent undersampling of histograms!
        if self.nzmax < len(self.z_samples) - 1:
            print("You're trying to integrate at lower resolution than supplied by the n(z) histograms. \n Increase nzmax>={:}! Aborting now...".format(len(self.z_samples) - 1))
            exit()
        # if that's the case, we want to integrate at histogram resolution and need to account for
        # the extra zero entry added
        elif self.nzmax == len(self.z_samples) - 1:
            self.nzmax = self.z_samples.shape[1]
            # requires that z-spacing is always the same for all bins...
            self.z_p = self.z_samples[0, :]
            print('Integrations performed at resolution of histogram! \n')
        # if we interpolate anyway at arbitrary resolution the extra 0 doesn't matter
        else:
            self.nzmax += 1
            self.z_p = np.linspace(self.z_samples.min(), self.z_samples.max(), self.nzmax)
            print('Integration performed at set nzmax={:} resolution! \n'.format(self.nzmax - 1))

        self.pz = np.zeros((self.nzmax, self.nzbins))
        self.pz_norm = np.zeros(self.nzbins, 'float64')
        self.splines_pz = []
        for zbin in xrange(self.nzbins):
                # we assume that the z-spacing is the same for each histogram
                spline_pz = itp.interp1d(self.z_samples[zbin, :], self.hist_samples[zbin, :], kind=self.type_redshift_interp)
                self.splines_pz.append(spline_pz)
                mask_min = self.z_p >= self.z_samples[zbin, :].min()
                mask_max = self.z_p <= self.z_samples[zbin, :].max()
                mask = mask_min & mask_max
                # points outside the z-range of the histograms are set to 0!
                #self.pz[mask, zbin] = itp.splev(self.z_p[mask], spline_pz)
                self.pz[mask, zbin] = spline_pz(self.z_p[mask])
                # Normalize selection functions
                dz = self.z_p[1:] - self.z_p[:-1]
                self.pz_norm[zbin] = np.sum(0.5 * (self.pz[1:, zbin] + self.pz[:-1, zbin]) * dz)

        self.zmax = self.z_p.max()
        self.need_cosmo_arguments(data, {'z_max_pk': self.zmax})

        # Initialize the BandPowers module from CosmoSIS:
        config_theory = {'bandpowers': {'theta_min': self.theta_min,
                                        'theta_max': self.theta_max,
                                        'l_min': self.ell_bin_min,
                                        'l_max': self.ell_bin_max,
                                        'nBands': self.nbins,
                                        'type': 'cosmic_shear_e',
                                        'Apodise': int(self.apodise),
                                        'Delta_x': self.delta_x,
                                        'Response_function_type': self.response_function,
                                        'Analytic': int(self.analytic),
                                        'input_section_name': 'shear_cl',
                                        'Output_FolderName': os.path.join(self.data_directory, 'BandPower_outputs')
                                        }}

        if hasattr(self, 'ell_bin_file'):
            config_theory.update({'l_min_l_max_file': self.ell_bin_file})

        # needs to point down to one of the '.so' files in '../kcap/cosebis/':
        self.theory_module = cosmosis.runtime.module.Module(module_name='bandpowers',
                                            file_path=os.path.join(self.kcap_directory,
                                            'cosebis/libbandpower.so'))

        self.theory_module.setup(dict_to_datablock(config_theory))

        return

    def loglkl(self, cosmo, data):

        theory_vec = self.cosmo_calculations(cosmo, data)
        
        if self.write_out_theory:
            fname = os.path.join(self.data_directory, self.theory_file)
            # for now we just dump the theory vector with no further details,
            # ASCII style...
            np.savetxt(fname, theory_vec)
            print('Saved theory vector to: \n {:} \n'.format(fname))
            print('Aborting run now. \n Set flag "write_out_theory = False" for likelihood evaluations \n and double-check your data vector!')
            exit()

        diff_vec = self.data_vec - theory_vec
        # apply mask:
        diff_vec = diff_vec

        # Don't invert that matrix!
        #chi2 = difference_vector.T.dot(inv_cov_sliced.dot(difference_vector))
        # this is for running smoothly with MultiNest
        # (in initial checking of prior space, there might occur weird solutions)
        if np.isinf(diff_vec).any() or np.isnan(diff_vec).any():
            chi2 = 2e12
        else:
            # don't invert that matrix...
            # use the Cholesky decomposition instead:
            yt = solve_triangular(self.cholesky_transform, diff_vec, lower=True)
            chi2 = yt.dot(yt)

        # enforce Gaussian priors on NUISANCE parameters if requested:
        if self.use_gaussian_prior_for_nuisance:

            for idx_nuisance, nuisance_name in enumerate(self.gaussian_prior_name):

                scale = data.mcmc_parameters[nuisance_name]['scale']
                chi2 += (data.mcmc_parameters[nuisance_name]['current'] * scale - self.gaussian_prior_center[idx_nuisance])**2 / self.gaussian_prior_sigma[idx_nuisance]**2

        #dt = timer() - t0
        #print('Time for one likelihood evaluation: {:.6f}s.'.format(dt))

        return -0.5 * chi2

    def load_data_file(self):
        """
        Function to load and process the data FITS file
        """

        # we really only need to load the n(z) data explicitly, the rest
        # is handled by the scale_cuts module
        data_tables = fits.open(os.path.join(self.data_directory, self.data_file))
        '''
        # this is only needed in the 2cosmos likelihood (for creating a mask):
        estimator = data_tables['PeeE'].data
        data_vec = np.asarray(estimator['VALUE'])
        covmat = np.asarray(data_tables['COVMAT'].data)
        '''

        # read number of NZ_SOURCE bins from HEADER
        # This does not seem to work and was a lucky hit...
        #self.nzbins = int(data_tables[2].header['N_ZBIN_1'])
        
        # define also number of unique z-bin correlations:
        self.nzcorrs = self.nzbins * (self.nzbins + 1) // 2

        nofz = data_tables['NZ_SOURCE'].data

        z_samples = []
        hist_samples = []
        for zbin in xrange(self.nzbins):
            zptemp= nofz['Z_MID']
            hist_pz = nofz['BIN{:}'.format(zbin + 1)]
            z_samples += [np.concatenate((np.zeros(1), zptemp))]
            hist_samples += [np.concatenate((np.zeros(1), hist_pz))]

        return np.asarray(z_samples), np.asarray(hist_samples)

    def baryon_feedback_bias_sqr(self, k, z, A_bary=1.):
        """

        Fitting formula for baryon feedback after equation 10 and Table 2 from J. Harnois-Deraps et al. 2014 (arXiv.1407.4301)

        """

        # k is expected in h/Mpc and is divided in log by this unit...
        x = np.log10(k)

        a = 1. / (1. + z)
        a_sqr = a * a

        constant = {'AGN':   {'A2': -0.11900, 'B2':  0.1300, 'C2':  0.6000, 'D2':  0.002110, 'E2': -2.0600,
                              'A1':  0.30800, 'B1': -0.6600, 'C1': -0.7600, 'D1': -0.002950, 'E1':  1.8400,
                              'A0':  0.15000, 'B0':  1.2200, 'C0':  1.3800, 'D0':  0.001300, 'E0':  3.5700},
                    'REF':   {'A2': -0.05880, 'B2': -0.2510, 'C2': -0.9340, 'D2': -0.004540, 'E2':  0.8580,
                              'A1':  0.07280, 'B1':  0.0381, 'C1':  1.0600, 'D1':  0.006520, 'E1': -1.7900,
                              'A0':  0.00972, 'B0':  1.1200, 'C0':  0.7500, 'D0': -0.000196, 'E0':  4.5400},
                    'DBLIM': {'A2': -0.29500, 'B2': -0.9890, 'C2': -0.0143, 'D2':  0.001990, 'E2': -0.8250,
                              'A1':  0.49000, 'B1':  0.6420, 'C1': -0.0594, 'D1': -0.002350, 'E1': -0.0611,
                              'A0': -0.01660, 'B0':  1.0500, 'C0':  1.3000, 'D0':  0.001200, 'E0':  4.4800}}

        A_z = constant[self.baryon_model]['A2']*a_sqr+constant[self.baryon_model]['A1']*a+constant[self.baryon_model]['A0']
        B_z = constant[self.baryon_model]['B2']*a_sqr+constant[self.baryon_model]['B1']*a+constant[self.baryon_model]['B0']
        C_z = constant[self.baryon_model]['C2']*a_sqr+constant[self.baryon_model]['C1']*a+constant[self.baryon_model]['C0']
        D_z = constant[self.baryon_model]['D2']*a_sqr+constant[self.baryon_model]['D1']*a+constant[self.baryon_model]['D0']
        E_z = constant[self.baryon_model]['E2']*a_sqr+constant[self.baryon_model]['E1']*a+constant[self.baryon_model]['E0']

        # only for debugging; tested and works!
        #print 'AGN: A2=-0.11900, B2= 0.1300, C2= 0.6000, D2= 0.002110, E2=-2.0600'
        #print self.baryon_model+': A2={:.5f}, B2={:.5f}, C2={:.5f}, D2={:.5f}, E2={:.5f}'.format(constant[self.baryon_model]['A2'], constant[self.baryon_model]['B2'], constant[self.baryon_model]['C2'],constant[self.baryon_model]['D2'], constant[self.baryon_model]['E2'])

        # original formula:
        #bias_sqr = 1.-A_z*np.exp((B_z-C_z)**3)+D_z*x*np.exp(E_z*x)
        # original formula with a free amplitude A_bary:
        bias_sqr = 1. - A_bary * (A_z * np.exp((B_z * x - C_z)**3) - D_z * x * np.exp(E_z * x))

        return bias_sqr

    def get_IA_factor(self, z, linear_growth_rate, rho_crit, Omega_m, small_h, amplitude, exponent):

        const = 5e-14 / small_h**2 # Mpc^3 / M_sol

        # arbitrary convention
        z0 = 0.3
        #print utils.growth_factor(z, self.Omega_m)
        #print self.rho_crit
        factor = -1. * amplitude * const * rho_crit * Omega_m / linear_growth_rate * ((1. + z) / (1. + z0))**exponent

        return factor

    def get_critical_density(self, small_h):
        """
        The critical density of the Universe at redshift 0.

        Returns
        -------
        rho_crit in solar masses per cubic Megaparsec.

        """

        # yay, constants...
        Mpc_cm = 3.08568025e24 # cm
        M_sun_g = 1.98892e33 # g
        G_const_Mpc_Msun_s = M_sun_g * (6.673e-8) / Mpc_cm**3.
        H100_s = 100. / (Mpc_cm * 1.0e-5) # s^-1

        rho_crit_0 = 3. * (small_h * H100_s)**2. / (8. * np.pi * G_const_Mpc_Msun_s)

        return rho_crit_0

    def get_matter_power_spectrum(self, r, z, cosmo, data):

        # Get power spectrum P(k=l/r,z(r)) from cosmological module
        pk = np.zeros((self.nells, self.nzmax), 'float64')
        pk_lin = np.zeros((self.nells, self.nzmax), 'float64')

        k_max_in_inv_Mpc = self.k_max_h_by_Mpc * cosmo.h()
        for idx_ell in xrange(self.nells):
            for idx_z in xrange(self.nzmax):
                # standard Limber approximation:
                #k = ells[idx_ell] / r[idx_z]
                # extended Limber approximation (cf. LoVerde & Afshordi 2008):
                k_in_inv_Mpc = (self.ells[idx_ell] + 0.5) / r[idx_z]
                if k_in_inv_Mpc > k_max_in_inv_Mpc:
                    pk_dm = 0.
                    pk_lin_dm = 0.
                else:
                    pk_dm = cosmo.pk(k_in_inv_Mpc, z[idx_z])
                    pk_lin_dm = cosmo.pk_lin(k_in_inv_Mpc, z[idx_z])

                if 'A_bary' in data.mcmc_parameters:
                    A_bary = data.mcmc_parameters['A_bary']['current'] * data.mcmc_parameters['A_bary']['scale']
                    pk[idx_ell, idx_z] = pk_dm * self.baryon_feedback_bias_sqr(k_in_inv_Mpc / cosmo.h(), z[idx_z], A_bary=A_bary)
                    # don't apply the baryon feedback model to the linear Pk!
                    #pk_lin[idx_ell, idx_z] = pk_lin_dm * self.baryon_feedback_bias_sqr(k_in_inv_Mpc / cosmo.h(), z[idx_z], A_bary=A_bary)
                else:
                    pk[idx_ell, idx_z] = pk_dm
                    pk_lin[idx_ell, idx_z] = pk_lin_dm

        return pk, pk_lin

    def get_lensing_kernel(self, r, pr):
        """
        Compute function g_i(r), that depends on r and the bin
        g_i(r) = 2r(1+z(r)) int_r^+\infty drs p_r(rs) (rs-r)/rs
        """

        g = np.zeros((self.nzmax, self.nzbins), 'float64')
        for Bin in xrange(self.nzbins):
            # shift only necessary if z[0] = 0
            for nr in xrange(1, self.nzmax - 1):
            #for nr in xrange(self.nzmax - 1):
                fun = pr[nr:, Bin] * (r[nr:] - r[nr]) / r[nr:]
                g[nr, Bin] = np.sum(0.5*(fun[1:] + fun[:-1]) * (r[nr+1:] - r[nr:-1]))
                g[nr, Bin] *= 2. * r[nr] * (1. + self.z_p[nr])

        return g

    def get_shear_power_spectrum(self, cosmo, data):
        """
        Function to calculate angular shear-shear power spectra, Cls.
        """

        # Omega_m contains all species!
        Omega_m = cosmo.Omega_m()
        small_h = cosmo.h()

        # needed for IA modelling:
        if ('A_IA' in data.mcmc_parameters) and ('exp_IA' in data.mcmc_parameters):
            amp_IA = data.mcmc_parameters['A_IA']['current'] * data.mcmc_parameters['A_IA']['scale']
            exp_IA = data.mcmc_parameters['exp_IA']['current'] * data.mcmc_parameters['exp_IA']['scale']
            intrinsic_alignment = True
        elif ('A_IA' in data.mcmc_parameters) and ('exp_IA' not in data.mcmc_parameters):
            amp_IA = data.mcmc_parameters['A_IA']['current'] * data.mcmc_parameters['A_IA']['scale']
            # redshift-scaling is turned off:
            exp_IA = 0.

            intrinsic_alignment = True
        else:
            intrinsic_alignment = False

        # One wants to obtain here the relation between z and r, this is done
        # by asking the cosmological module with the function z_of_r
        r, dzdr = cosmo.z_of_r(self.z_p)

        # Compute now the selection function p(r) = p(z) dz/dr normalized
        # to one. The np.newaxis helps to broadcast the one-dimensional array
        # dzdr to the proper shape. Note that p_norm is also broadcasted as
        # an array of the same shape as p_z
        if (self.shift_n_z_by_D_z.any()):

            # correlate D_z shifts:
            D_z = np.zeros(self.nzbins)
            for zbin in xrange(self.nzbins):

                param_name = 'D_z{:}'.format(zbin + 1)
                if param_name in data.mcmc_parameters:
                    D_z[zbin] = data.mcmc_parameters[param_name]['current'] * data.mcmc_parameters[param_name]['scale']

            D_z_corr = self.L_matrix_D_z.dot(D_z)

            pz = np.zeros((self.nzmax, self.nzbins), 'float64')
            pz_norm = np.zeros(self.nzbins, 'float64')
            for zbin in xrange(self.nzbins):

                '''
                param_name = 'D_z{:}'.format(zbin + 1)
                if param_name in data.mcmc_parameters:
                    z_mod = self.z_p + data.mcmc_parameters[param_name]['current'] * data.mcmc_parameters[param_name]['scale']
                else:
                    z_mod = self.z_p
                '''
                z_mod = self.z_p + D_z_corr[zbin]
                spline_pz = self.splines_pz[zbin]
                mask_min = z_mod >= self.z_samples[zbin, :].min()
                mask_max = z_mod <= self.z_samples[zbin, :].max()
                mask = mask_min & mask_max
                # points outside the z-range of the histograms are set to 0!
                #pz[mask, zbin] = itp.splev(z_mod[mask], spline_pz)
                pz[mask, zbin] = spline_pz(z_mod[mask])
                # Normalize selection functions
                dz = self.z_p[1:] - self.z_p[:-1]
                pz_norm[zbin] = np.sum(0.5 * (pz[1:, zbin] + pz[:-1, zbin]) * dz)

            pr = pz * (dzdr[:, np.newaxis] / pz_norm)

        else:
            # use fiducial dn/dz loaded in the __init__:
            pr = self.pz * (dzdr[:, np.newaxis] / self.pz_norm)

        # get linear growth rate if IA are modelled:
        if intrinsic_alignment:
            rho_crit = self.get_critical_density(small_h)
            # derive the linear growth factor D(z)
            linear_growth_rate = np.zeros_like(self.z_p)
            #print self.redshifts
            for idx_z, z in enumerate(self.z_p):
                try:
                    # for CLASS ver >= 2.6:
                    linear_growth_rate[idx_z] = cosmo.scale_independent_growth_factor(z)
                except:
                    # my own function from private CLASS modification:
                    linear_growth_rate[idx_z] = cosmo.growth_factor_at_z(z)
            # normalize to unity at z=0:
            try:
                # for CLASS ver >= 2.6:
                linear_growth_rate /= cosmo.scale_independent_growth_factor(0.)
            except:
                # my own function from private CLASS modification:
                linear_growth_rate /= cosmo.growth_factor_at_z(0.)

        g = self.get_lensing_kernel(r, pr)
        pk, pk_lin = self.get_matter_power_spectrum(r, self.z_p, cosmo, data)

        Cl_integrand = np.zeros((self.nzmax, self.nzcorrs), 'float64')
        Cl = np.zeros((self.nzcorrs, self.nells), 'float64')

        Cl_GG_integrand = np.zeros_like(Cl_integrand)
        Cl_GG = np.zeros_like(Cl)

        if intrinsic_alignment:
            Cl_II_integrand = np.zeros_like(Cl_integrand)
            Cl_II = np.zeros_like(Cl)

            Cl_GI_integrand = np.zeros_like(Cl_integrand)
            Cl_GI = np.zeros_like(Cl)

        list_cl_keys = []
        dr = r[1:] - r[:-1]
        # Start loop over l for computation of C_l^shear
        for il in xrange(self.nells):
            # find Cl_integrand = (g(r) / r)**2 * P(l/r,z(r))
            for Bin1 in xrange(self.nzbins):
                for Bin2 in xrange(Bin1, self.nzbins):
                    if il == 0:
                        list_cl_keys += ['bin_{:}_{:}'.format(Bin2 + 1, Bin1 + 1)]
                    Cl_GG_integrand[1:, self.__one_dim_index(Bin1,Bin2)] = g[1:, Bin1] * g[1:, Bin2] / r[1:]**2 * pk[il, 1:]
                    if intrinsic_alignment:
                        factor_IA = self.get_IA_factor(self.z_p, linear_growth_rate, rho_crit, Omega_m, small_h, amp_IA, exp_IA) #/ self.dzdr[1:]

                        if self.use_linear_pk_for_IA:
                            # this term (II) uses the linear matter power spectrum P_lin(k, z)
                            Cl_II_integrand[1:, self.__one_dim_index(Bin1, Bin2)] = pr[1:, Bin1] * pr[1:, Bin2] * factor_IA[1:]**2 / r[1:]**2 * pk_lin[il, 1:]
                            # this term (GI) uses sqrt(P_lin(k, z) * P_nl(k, z))
                            Cl_GI_integrand[1:, self.__one_dim_index(Bin1, Bin2)] = (g[1:, Bin1] * pr[1:, Bin2] + g[1:, Bin2] * pr[1:, Bin1]) * factor_IA[1:] / r[1:]**2 * np.sqrt(pk_lin[il, 1:] * pk[il, 1:])
                        else:
                            # both II and GI terms use the non-linear matter power spectrum P_nl(k, z)
                            Cl_II_integrand[1:, self.__one_dim_index(Bin1, Bin2)] = pr[1:, Bin1] * pr[1:, Bin2] * factor_IA[1:]**2 / r[1:]**2 * pk[il, 1:]
                            Cl_GI_integrand[1:, self.__one_dim_index(Bin1, Bin2)] = (g[1:, Bin1] * pr[1:, Bin2] + g[1:, Bin2] * pr[1:, Bin1]) * factor_IA[1:] / r[1:]**2 * pk[il, 1:]

            # Integrate over r to get C_l^shear_ij = P_ij(l)
            # C_l^shear_ij = 9/16 Omega0_m^2 H_0^4 \sum_0^rmax dr (g_i(r)
            # g_j(r) /r**2) P(k=l/r,z(r)) dr
            # It is then multiplied by 9/16*Omega_m**2
            # and then by (h/2997.9)**4 to be dimensionless
            # (since P(k)*dr is in units of Mpc**4)
            for Bin in xrange(self.nzcorrs):
                Cl_GG[Bin, il] = np.sum(0.5 * (Cl_GG_integrand[1:, Bin] + Cl_GG_integrand[:-1, Bin]) * dr)
                Cl_GG[Bin, il] *= 9. / 16. * Omega_m**2
                Cl_GG[Bin, il] *= (small_h / 2997.9)**4

                if intrinsic_alignment:
                    Cl_II[Bin, il] = np.sum(0.5 * (Cl_II_integrand[1:, Bin] + Cl_II_integrand[:-1, Bin]) * dr)

                    Cl_GI[Bin, il] = np.sum(0.5 * (Cl_GI_integrand[1:, Bin] + Cl_GI_integrand[:-1, Bin]) * dr)
                    # here we divide by 4, because we get a 2 from g(r)!
                    Cl_GI[Bin, il] *= 3. / 4. * Omega_m
                    Cl_GI[Bin, il] *= (small_h / 2997.9)**2

        if intrinsic_alignment:
            Cl = Cl_GG + Cl_GI + Cl_II
        else:
            Cl = Cl_GG

        if self.write_out_Cls:
            #Cls_out = np.zeros((self.nzcorrs + 1, self.nells), 'float64')
            Cls_out = self.ells
            fname = os.path.join(self.data_directory, 'Cls_tot.txt')
            header = 'ells, '
            for idx in xrange(self.nzcorrs):
                header += list_cl_keys[idx] + ', '
                Cls_out = np.column_stack((Cls_out, Cl_GG[idx, :]))
            header = header[:-2]
            np.savetxt(fname, Cls_out, header=header)
            print('Saved Cls to: \n {:} \n'.format(fname))

        return Cl, list_cl_keys

    def get_theory_vec(self, Cls, Cl_keys, data):

        # create input dict for datablock:
        input_theory = {}
        # we need to set the number of nzbins:
        # just setting 'nbin' doesn't seem to work...
        #input_cosebis['shear_cl'] = {'nbin': 5}
        # this produces some output at last...
        input_theory['shear_cl'] = {'nbin_a': self.nzbins, 'nbin_b': self.nzbins}
        # now add the vals for 'ell' and 'bin_1_1', 'bin_2_1', ... 'bin_n_n'
        input_theory['shear_cl'].update({'ell': self.ells})
        input_theory['shear_cl'].update(dict(zip(Cl_keys, Cls[:, :])))

        datablock = dict_to_datablock(input_theory)
        #print(block_cosebis.keys())
        self.theory_module.execute(datablock)
        
        # silence the scale_cuts module during likelihood evaluations:
        block_print()
        self.scale_cuts_module.execute(datablock)
        # re-enable print-statements again:
        enable_print()
        
        theory_vec = np.asarray(datablock['likelihood_bp', 'theory'])

        return theory_vec

    def cosmo_calculations(self, cosmo, data):

        Cls, Cl_keys = self.get_shear_power_spectrum(cosmo, data)
        theory_vec = self.get_theory_vec(Cls, Cl_keys, data)

        return theory_vec

    def __one_dim_index(self, Bin1, Bin2):
        """
        This function is used to convert 2D sums over the two indices (Bin1, Bin2)
        of an N*N symmetric matrix into 1D sums over one index with N(N+1)/2
        possible values.
        """

        if Bin1 <= Bin2:
            return Bin2 + self.nzbins * Bin1 - (Bin1 * (Bin1 + 1)) // 2
        else:
            return Bin1 + self.nzbins * Bin2 - (Bin2 * (Bin2 + 1)) // 2