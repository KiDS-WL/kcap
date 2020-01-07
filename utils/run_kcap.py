import os
import shutil
import warnings
import configparser
import argparse

import numpy as np

import sys
sys.path.append("kcap")
import cosmosis_utils

from cosmosis.datablock import names


# KCAP_PATH = os.path.join(os.path.dirname(__file__), "..")
KCAP_PATH = "."
BOSS_PATH = os.path.join(KCAP_PATH, "../kcap_boss_module")
CSL_PATH = os.path.join(KCAP_PATH, "cosmosis-standard-library")
COSEBIS_PATH = os.path.join(KCAP_PATH, "cosebis")
COSEBIS_OUTPUTS = COSEBIS_PATH

SCALE_CUT_PATH = os.path.join(KCAP_PATH, "modules/scale_cuts")
MOCK_DIR = os.path.join(KCAP_PATH, "data/gaussian_mocks/KV450/")

class K1000Pipeline:
    def __init__(self, **kwargs):

        self.full_config = self.create_config(**kwargs)
        self.full_values, self.full_priors = self.create_values()

        self.pipelines = {"EE_nE_w" :         {"cut_modules" : ["cl2xi_shear", "cl2xi_ggl", "bin_xi_plus", "bin_xi_minus", "bin_xi_ggl"],
                                               "cut_keys"    : [],},
                          "EE_nE_w_mocks" :   {"cut_modules" : ["correlated_dz_priors", "source_photoz_bias",
                                                                "cl2xi_shear", "cl2xi_ggl", "bin_xi_plus", "bin_xi_minus", "bin_xi_ggl"],
                                               "cut_values"  : ["nofz_shifts"],
                                               "sample"      : False,},
                          
                          "EE" :              {"cut_modules" : ["wedges", "BOSS_like", "approximate_P_gm",
                                                                "load_lens_nz",
                                                                "cl2xi_shear", "cl2xi_ggl", "bin_xi_plus", "bin_xi_minus", "bin_xi_ggl",
                                                                "bandpower_ggl"],
                                               "cut_keys"    : [("projection", "position-shear")],
                                               "set_keys"    : [("scale_cuts", "use_stats", "PeeE")],
                                               "cut_values"  : ["bias_parameters"]},
                          "EE_mocks" :        {"cut_modules" : ["correlated_dz_priors", "source_photoz_bias",
                                                                "wedges", "BOSS_like", "approximate_P_gm",
                                                                "load_lens_nz",
                                                                "cl2xi_shear", "cl2xi_ggl", "bin_xi_plus", "bin_xi_minus", "bin_xi_ggl",
                                                                "bandpower_ggl"],
                                               "cut_keys"    : [("projection", "position-shear")],
                                               "set_keys"    : [("scale_cuts", "use_stats", "PeeE")],
                                               "cut_values"  : ["nofz_shifts", "bias_parameters"],
                                               "sample"      : False,},

                          "EE_nE" :           {"cut_modules" : ["approximate_P_gm",
                                                                "cl2xi_shear", "cl2xi_ggl", "bin_xi_plus", "bin_xi_minus", "bin_xi_ggl",],},
                          "EE_nE_mocks" :     {"cut_modules" : ["correlated_dz_priors", "source_photoz_bias",
                                                                "wedges", "BOSS_like",
                                                                "cl2xi_shear", "cl2xi_ggl", "bin_xi_plus", "bin_xi_minus", "bin_xi_ggl",],
                                               "cut_values"  : ["nofz_shifts"],
                                               "sample"      : False,},
                         }

    def choose_pipeline(self, name=None, pipeline=None):
        if name is not None:
            pipeline = self.pipelines[name]

        config = {**self.full_config}
        values = {**self.full_values}
        priors = {**self.full_priors}

        if "cut_modules" in pipeline:
            for mod in pipeline["cut_modules"]:
                del config[mod]

        if "cut_keys" in pipeline:
            for mod, key in pipeline["cut_keys"]:
                del config[mod][key]

        if "set_keys" in pipeline:
            for mod, key, val in pipeline["set_keys"]:
                config[mod][key] = val

        if "cut_values" in pipeline:
            for val in pipeline["cut_values"]:
                del values[val]
                priors.pop(val, None)

        if not pipeline.get("sample", True):
            # Not sampling, i.e., create mocks
            for section in values.keys():
                for key, value in values[section].items():
                    if isinstance(value, (list, tuple)):
                        # Take starting point of range
                        values[section][key] = value[1]
            priors = {}


        self.config = config
        self.values = values
        self.priors = priors

    def run_pipeline(self, params=None):
        if params is None:
            params = self.values

        pipeline = cosmosis_utils.create_pipeline(self.config)

        block = cosmosis_utils.dict_to_datablock(params)
        pipeline(block)

        for k, v in block.keys():
            print(k, v)

        return block

    def build_ini(self):
        ini = configparser.ConfigParser()
        ini.read_dict(self.flatten_config({**self.sampling_config, **self.config}))

        values_ini = configparser.ConfigParser()
        values_ini.read_dict(self.flatten_config(self.values))

        priors_ini = configparser.ConfigParser()
        priors_ini.read_dict(self.flatten_config(self.priors))

        print(cosmosis_utils.config_to_string(ini))
        print(cosmosis_utils.config_to_string(values_ini))
        print(cosmosis_utils.config_to_string(priors_ini))

        return ini, values_ini, priors_ini

    def prepare_run_directory(self, output_root_dir, run_name, 
                                    ini, values_ini, priors_ini,
                                    verbose=True):
        run_dir = os.path.join(output_root_dir, run_name, "config")

        if verbose: print(f"Creating run directory for {run_name}: {run_dir}")
        os.makedirs(run_dir)

        # Write ini file
        ini.write(os.path.join(run_dir, f"{run_name}.ini"))
        
        # Write values file
        values_ini.write(os.path.join(run_dir, f"{run_name}_values.ini"))

        # Write priors file
        if len(priors_ini) > 0:
            priors_ini.write(os.path.join(run_dir, f"{run_name}_priors.ini"))
        
        # Write data file

        # Write any other file used.
        # log git commit
        # Create tarball script.

    def create_config(self, KiDS_data_cov_file,
                            BOSS_data_files,
                            BOSS_cov_files,
                            KiDS_mock_output_file,
                            BOSS_mock_output_files,
                            source_nz_files,
                            lens_nz_files,
                            dz_covariance_file,
                            source_nz_sample,
                            lens_nz_sample,
                            BOSS_like_name,
                            bands_range,
                            points_range,
                            num_ell,
                            z_eff,
                            xi_theta_min,
                            xi_theta_max,
                            xi_n_theta,
                            bandpower_ell_min,
                            bandpower_ell_max,
                            bandpower_n_bin,
                            bandpower_theta_min,
                            bandpower_theta_max,
                            bandpower_apodise,
                            bandpower_Delta_x,
                            used_stats,
                            halofit_version,
                            create_mocks,
                            noisy_mocks):
        
        config = {  "sample_ln_As"       : {"file" : os.path.join(KCAP_PATH,
                                                            "utils/sample_ln_As.py"),},

                    "correlated_dz_priors":{"file" : os.path.join(KCAP_PATH,
                                                            "utils/correlated_priors.py"),
                                            "uncorrelated_parameters" : " ".join(["nofz_shifts/p_1",
                                                                                "nofz_shifts/p_2",
                                                                                "nofz_shifts/p_3",
                                                                                "nofz_shifts/p_4",
                                                                                "nofz_shifts/p_5",]),
                                            "output_parameters"       : " ".join(["nofz_shifts/bias_1",
                                                                                "nofz_shifts/bias_2",
                                                                                "nofz_shifts/bias_3",
                                                                                "nofz_shifts/bias_4",
                                                                                "nofz_shifts/bias_5",]),
                                            "covariance"              : dz_covariance_file},
                    
                    "one_parameter_hmcode":{"file" : os.path.join(KCAP_PATH,
                                                            "utils/one_parameter_hmcode.py"),
                                            "a_0" : 0.98,
                                            "a_1" : -0.12},

                    "camb"               : {"file" : os.path.join(CSL_PATH, 
                                                                "boltzmann/pycamb/camb_interface.py"),
                                            "do_reionization"    : False,
                                            "mode"               : "transfer",
                                            "nonlinear"          : "pk",
                                            "halofit_version"    : halofit_version,
                                            "neutrino_hierarchy" : "normal",
                                            "kmax"               : 20.0,
                                            "zmid"               : 2.0,
                                            "nz_mid"             : 100,
                                            "zmax"               : 6.0,
                                            "nz"                 : 150,
                                            "background_zmax"    : 6.0,
                                            "background_zmin"    : 0.0,
                                            "background_nz"      : 6000,
                                            },

                    "wedges"     :         {"file" : os.path.join(BOSS_PATH, 
                                                                "python_interface/cosmosis_module.py"),
                                            "window_file" : os.path.join(BOSS_PATH,
                                                                        "CosmoMC_BOSS/data/BOSS.DR12_windows.txt"),
                                            "bands_file" : os.path.join(BOSS_PATH,
                                                                        "CosmoMC_BOSS/data/BOSS.DR12_rbands.txt"),
                                            "bands_range"         : bands_range,
                                            "points_range"        : points_range,
                                            "num_ell"             : num_ell,
                                            "z_eff"               : z_eff,
                                            "output_section_name" : "xi_wedges",
                                            "use_growth"          : False,
                                            "local_lag_g2"        : True,
                                            "local_lag_g3"        : False,},

                    "approximate_P_gm" :   {"file" : os.path.join(KCAP_PATH,
                                                                "modules/P_gm_approx/p_gm_approx_interface.py"),
                                            "b2_coefficient_file" : os.path.join(KCAP_PATH,
                                                                        "modules/P_gm_approx/parameter_gridfit_b2.dat"),
                                            "g2_coefficient_file" : os.path.join(KCAP_PATH,
                                                                        "modules/P_gm_approx/parameter_gridfit_g2.dat"),
                                            "g3_coefficient_file" : os.path.join(KCAP_PATH,
                                                                        "modules/P_gm_approx/parameter_gridfit_g3.dat"),
                                            "z_sep"               : 0.5,},

                    "extrapolate_power" :  {"file" : os.path.join(CSL_PATH, 
                                                        "boltzmann/extrapolate/extrapolate_power.py"),
                                            "kmax" : 500.0},

                    "load_source_nz" :     {"file" : os.path.join(CSL_PATH, 
                                                            "number_density/load_nz/load_nz.py"),
                                            "filepath" : " ".join(source_nz_files),
                                            "histogram" : True,
                                            "output_section" : f"nz_{source_nz_sample}"},
                        
                    "load_lens_nz"       : {"file" : os.path.join(CSL_PATH, 
                                                            "number_density/load_nz/load_nz.py"),
                                            "filepath" : " ".join(lens_nz_files),
                                            "histogram" : True,
                                            "output_section" : f"nz_{lens_nz_sample}"},

                    "source_photoz_bias" : {"file" : os.path.join(CSL_PATH, 
                                                        "number_density/photoz_bias/photoz_bias.py"),
                                            "mode" : "additive",
                                            "sample" : f"nz_{source_nz_sample}",
                                            "bias_section" : "nofz_shifts",
                                            "interpolation" : "linear"},

                    "linear_alignment" :   {"file" : os.path.join(CSL_PATH, 
                                                        "intrinsic_alignments/la_model/linear_alignments_interface.py"),
                                            "method" : "bk_corrected"},

                    "projection" :         {"file" : os.path.join(CSL_PATH, 
                                                        "structure/projection/project_2d.py"),
                                            "ell_min" : 10.0,
                                            "ell_max" : 3.0e4,
                                            "n_ell"  : 400,
                                            "shear-shear" : f"{source_nz_sample}-{source_nz_sample}",
                                            "position-shear" : f"{lens_nz_sample}-{source_nz_sample}",
                                            "shear-intrinsic" : f"{source_nz_sample}-{source_nz_sample}",
                                            #"position-intrinsic" : f"{lens_nz_sample}-{source_nz_sample}",
                                            "intrinsic-intrinsic" : f"{source_nz_sample}-{source_nz_sample}",
                                            "verbose" : False,
                                            "get_kernel_peaks" : False},


                    "add_intrinsic" :      {"file" : os.path.join(CSL_PATH, 
                                                        "shear/add_intrinsic/add_intrinsic.py"),
                                            "position-shear" : False},

                    "cl2xi_shear" :        {"file" : os.path.join(CSL_PATH, 
                                                        "shear/cl_to_xi_nicaea/nicaea_interface.so"),
                                            "corr_type" : 0},
                    "cl2xi_ggl" :          {"file" : os.path.join(CSL_PATH, 
                                                        "shear/cl_to_xi_nicaea/nicaea_interface.so"),
                                            "corr_type" : 2},

                    "bin_xi_plus" :        {"file" : os.path.join(COSEBIS_PATH, "libxipm_binned.so"),
                                            "output_section_name" : "shear_xi_plus_binned", # (optional) the DEFAULT is xi_binned
                                            "input_section_name"  : "shear_xi_plus", # (optional) the DEFAULT depends on type
                                            "type"                : "plus", # please specify this otherwise as plus or minus DEFAULT is plus

                                            "theta_min"           : xi_theta_min,
                                            "theta_max"           : xi_theta_min,
                                            "nTheta"              : xi_n_theta,

                                            #"theta_min_max_filename" : "theta_bin_edges_file.ascii", # (optional) these are the edges of the theta plus bins,

                                            # InputNpair = InputNpair; (optional) a file containing the number of npair per finely binned thetas.
                                            # InputNpair_suffix = .ascii ; (optional) DEFAULT is empty
                                            # Column_theta = 0 ; (optional) which column in the file is theta? DEFAULT is 0
                                            # Column_Npair = 7 ; which column in the file is npair? DEFAULT is 7
                                            # nBins_in = 5 ; number of redshift bins, this needs to be given, otherwise will set weighted binning to just theta
                                            },
                    "bin_xi_minus" :       {"file" : os.path.join(COSEBIS_PATH, "libxipm_binned.so"),
                                            "output_section_name" : "shear_xi_minus_binned", # (optional) the DEFAULT is xi_binned
                                            "input_section_name"  : "shear_xi_minus", # (optional) the DEFAULT depends on type
                                            "type"                : "minus", # please specify this otherwise as plus or minus DEFAULT is plus

                                            "theta_min"           : xi_theta_min,
                                            "theta_max"           : xi_theta_min,
                                            "nTheta"              : xi_n_theta,},

                    "bin_xi_ggl" :         {"file" : os.path.join(COSEBIS_PATH, "libxipm_binned.so"),
                                            "output_section_name" : "galaxy_shear_xi_binned", # (optional) the DEFAULT is xi_binned
                                            "input_section_name"  : "galaxy_shear_xi", # (optional) the DEFAULT depends on type

                                            "theta_min"           : xi_theta_min,
                                            "theta_max"           : xi_theta_min,
                                            "nTheta"              : xi_n_theta,},

                    "bandpower_shear_e" :  {"file" : os.path.join(COSEBIS_PATH, "libbandpower.so"),
                                            "type"                   : "cosmic_shear_e",
                                            "Response_function_type" : "tophat",
                                            "Analytic"               : 1, #use analytic solution for g
                                            "output_section_name"    :  "bandpower_shear_e", # the DEFAULT is band_power
                                            #input_section_name = shear_cl ; the DEFAULT is shear_cl
                                            #l_min_max_file = l_min_max_file.ascii; a file with minimum and maximum values for ell. If it doesn't exist 
                                            # will look for l_min, l_max and nBands then do log binning between l_min and l_max.
                                            # if the file exists we will ignore l_min,l_max and nBands
                                            "l_min"                  : bandpower_ell_min,
                                            "l_max"                  : bandpower_ell_max,
                                            "nBands"                 : bandpower_n_bin,

                                            "Apodise"                : bandpower_apodise,
                                            "Delta_x"                : bandpower_Delta_x, # apodisation length in arcmins

                                            "theta_min"              : bandpower_theta_min,
                                            "theta_max"              : bandpower_theta_max,

                                            "Output_FolderName"      : os.path.join(COSEBIS_OUTPUTS, "BandPower_outputs/"),
                                            },
                    "bandpower_ggl" :      {"file" : os.path.join(COSEBIS_PATH, "libbandpower.so"),
                                            "type"                   : "ggl",
                                            "Response_function_type" : "tophat",
                                            "Analytic"               : 1, #use analytic solution for g
                                            "output_section_name"    :  "bandpower_galaxy_shear", # the DEFAULT is band_power
                                            #"input_section_name" = shear_cl ; the DEFAULT is shear_cl
                                            #l_min_max_file = l_min_max_file.ascii; a file with minimum and maximum values for ell. If it doesn't exist 
                                            # will look for l_min, l_max and nBands then do log binning between l_min and l_max.
                                            # if the file exists we will ignore l_min,l_max and nBands
                                            "l_min"                  : bandpower_ell_min,
                                            "l_max"                  : bandpower_ell_max,
                                            "nBands"                 : bandpower_n_bin,

                                            "Apodise"                : bandpower_apodise,
                                            "Delta_x"                : bandpower_Delta_x, # apodisation length in arcmins

                                            "theta_min"              : bandpower_theta_min,
                                            "theta_max"              : bandpower_theta_max,

                                            "Output_FolderName"      : os.path.join(COSEBIS_OUTPUTS, "BandPower_outputs/"),
                                            },

                    "scale_cuts"  :         {"file" : os.path.join(SCALE_CUT_PATH, "scale_cuts.py"),
                                            "output_section_name" : "theory_data_covariance",
                                            "data_and_covariance_fits_filename" : KiDS_data_cov_file,
                                            # Define the statistics to use & the scale cuts in the following two lines:
                                            #"scale_cuts_filename" : os.path.join(SCALE_CUT_PATH, "scale_cuts.ini"),
                                            #"scale_cuts_option" : "scale_cuts_none",
                                            "use_stats" : " ".join(used_stats),
    
                                            # Section names for data
                                            "xi_plus_extension_name" : "xiP",
                                            "xi_minus_extension_name" : "xiM",
                                            "bandpower_ggl_extension_name" : "PneE",
                                            "bandpower_e_cosmic_shear_extension_name" : "PeeE",
                                            "cosebis_extension_name" : "En",
    
                                            # Section names for theory
                                            "xi_plus_section_name" : "shear_xi_plus_binned",
                                            "xi_minus_section_name" : "shear_xi_minus_binned",
                                            "bandpower_ggl_section_name" : "bandpower_galaxy_shear",
                                            "bandpower_e_cosmic_shear_section_name" : "bandpower_shear_e",
                                            #cosebis_section_name = cosebis

                                            "simulate"  : create_mocks,
                                            "simulate_with_noise" : noisy_mocks,
                                            "mock_filename" : KiDS_mock_output_file,
                                            },


                    "BOSS_like"  : {"file" :  os.path.join(KCAP_PATH,
                                                        "utils/mini_BOSS_like.py"),
                                    "data_vector_file" : " ".join(BOSS_data_files),
                                    "covariance_file" : " ".join(BOSS_cov_files),
                                    "points_range"       : points_range,
                                    "like_name"          : BOSS_like_name,
                                    "keep_theory_vector" : True,},
                                    
                    "2x2pt_like" : {"file" : os.path.join(KCAP_PATH, "utils/mini_like.py"),
                                    "input_section_name" : "theory_data_covariance",
                                    "like_name"          : "2x2pt_like"},
        }
        return config

    def set_sampling_config(self, **kwargs):
        if not "likelihoods" in kwargs:
            kwargs["likelihoods"] = [m["like_name"] for m in self.config.values() if "like_name" in m]
        self.sampling_config = self.create_sampling_config(modules=self.config.keys(), **kwargs)

    def create_sampling_config(self, modules,
                                     likelihoods,
                                     derived_parameters,
                                     parameter_file,
                                     prior_file,
                                     verbose,
                                     debug,
                                     sampler_name,
                                     output_dir,
                                     run_name,
                                     max_iterations,
                                     resume=True,
                                     live_points=250,
                                     nested_sampling_tolerance=0.1,
                                     multinest_efficiency=0.8,
                                     multinest_const_efficiency=False,
                                     emcee_walker=80,
                                     emcee_covariance_file=None,
                                     ):
        config = {  "pipeline" :   {"modules"           : " ".join(modules),
                                    "values"            : parameter_file,
                                    "priors"            : prior_file,
                                    "likelihoods"       : " ".join(likelihoods),
                                    "extra_output"      : " ".join(derived_parameters),
                                    "quiet"             : "F" if verbose else "T",
                                    "timing"            : "T",
                                    "debug"             : "T" if debug else "F"},
            
                    "runtime"  :   {"sampler"           : sampler_name},
                  
                    "output"   :   {"filename"          : os.path.join(output_dir, f"samples_{run_name}.txt"),
                                    "format"            : "text"},
                  }

        samplers = {"multinest" :  {"max_iterations"    : max_iterations,
                                    "multinest_outfile_root" : os.path.join(output_dir, "multinest", f"multinest_{run_name}_"),
                                    "resume"            : "T" if resume else "F",
                                    "live_points"       : live_points,
                                    "efficiency"        : multinest_efficiency,
                                    "tolerance"         : nested_sampling_tolerance,
                                    "constant_efficiency" : "T" if multinest_const_efficiency else "F"},

                    "emcee"      : {"walkers"           : emcee_walker,
                                    "samples"           :  max_iterations,
                                    "covmat"            : emcee_covariance_file,
                                    "nsteps"            : 5},

                    "metropolis" : {"samples"           : max_iterations,
                                    "nsteps"            : 1},
                    }

        config[sampler_name] = samplers[sampler_name]

        return config

    def create_values(self):
        values = {"cosmological_parameters" :      {"omch2"       : [ 0.051,  0.13,      0.2],
                                                    "ombh2"       : [ 0.019,  0.0225,    0.026],
                                                    "h0"          : [ 0.64,   0.7,       0.82],
                                                    "n_s"         : [ 0.84,   0.97,      1.1],
                                                    "ln_1e10_A_s" : [ 1.5,    2.72,      4.0],
                                                    "omega_k"     :           0.0,
                                                    "w"           :          -1.0,
                                                    "mnu"         :           0.06,             #normal hierarchy
                                                    }, 
                "halo_model_parameters" :          {"A"           : [ 2.0,    2.6,       3.13]},

                "intrinsic_alignment_parameters" : {"A"           : [-6.0,    0.8,       6.0]},

                "nofz_shifts"                    : {"p_1"         : [-5.0,    0.0,       5.0],
                                                    "p_2"         : [-5.0,    0.0,       5.0],
                                                    "p_3"         : [-5.0,    0.0,       5.0],
                                                    "p_4"         : [-5.0,    0.0,       5.0],
                                                    "p_5"         : [-5.0,    0.0,       5.0],},

                "bias_parameters"         :        {"b1_bin_1"    : [0.5,     2.1,       9.0],
                                                    "b2_bin_1"    : [-4.0,    0.2,       8.0],
                                                    "gamma3_bin_1": [-8.0,    0.9,       8.0],
                                                    "a_vir_bin_1" : [0.0,     3.8,      12.0],

                                                    "b1_bin_2"    : [0.5,     2.3,       9.0],
                                                    "b2_bin_2"    : [-1.0,    0.5,       8.0],
                                                    "gamma3_bin_2": [-8.0,    0.1,       8.0],
                                                    "a_vir_bin_2" : [0.0,     3.0,      12.0],}}

        priors = {"nofz_shifts"                  : {"p_1"         : ["gaussian",    0.0,       1.0],
                                                    "p_2"         : ["gaussian",    0.0,       1.0],
                                                    "p_3"         : ["gaussian",    0.0,       1.0],
                                                    "p_4"         : ["gaussian",    0.0,       1.0],
                                                    "p_5"         : ["gaussian",    0.0,       1.0],}}

        return values, priors
        
    @staticmethod
    def flatten_config(values):
        values = {**values}
        for section, d in values.items():
            for key, value in d.items():
                if isinstance(value, (list, tuple)):
                    d[key] = " ".join([str(v) for v in value])
                elif value is True:
                    d[key] = "T"
                elif value is False:
                    d[key] = "F"
        return values
    
    
    def create_mock_BOSS(self, block, BOSS_like_name, num_ell, noisy_mocks=True):
        mock_data_BOSS = {}
        for b in range(len(z_eff)):
            if noisy_mocks:
                d = block[names.data_vector, BOSS_like_name + f"_simulation_bin_{b+1}"].reshape(num_ell, -1)
            else:
                d = block[names.data_vector, BOSS_like_name + f"_theory_bin_{b+1}"].reshape(num_ell, -1)
            err = np.sqrt(np.diag(block[names.data_vector, BOSS_like_name + f"_cov_bin_{b+1}"])).reshape(num_ell, -1)
            assert s.shape == d[0].shape
            mock_data_BOSS[b+1] = {"s"         : s,
                                f"w1"       : d[0],
                                f"w2"       : d[1],
                                f"w3"       : d[2],
                                f"w1_err"   : err[0],
                                f"w2_err"   : err[1],
                                f"w3_err"   : err[2],}
        return mock_data_BOSS


def create_BOSS_data_file(mock_data, filename):
    np.savetxt(filename, np.vstack((mock_data["s"], 
                                    mock_data["w1"], mock_data["w1_err"],
                                    mock_data["w2"], mock_data["w2_err"],
                                    mock_data["w3"], mock_data["w3_err"])).T,
                comments="#   s/(Mpc/h)          xi_3w,1       sigma xi_3w,1         xi_3w,2        sigma xi_3w,2        xi_3w,3        sigma xi_3w,3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--run-type", required=True)

    parser.add_argument("--mock-dir")
    parser.add_argument("--noiseless-mocks", action="store_true")

    parser.add_argument("--run-dir")
    parser.add_argument("--run-name")

    parser.add_argument("--halofit-version", default="mead")

    parser.add_argument("--sampler", default="multinest")

    args = parser.parse_args()

    pipeline_name = args.run_type
    create_mocks = args.mock_dir is not None
    noisy_mocks = not args.noiseless_mocks

    if create_mocks:
        mock_dir = args.mock_dir
    else:
        run_dir = args.run_dir
        run_name = args.run_name
        mock_dir = os.path.join(run_dir, run_name, "data")
        if not run_dir or not run_name:
            parser.error("If not creating mocks, --run-dir and --run-name need to be specified.")

    KiDS_data_cov_file = os.path.join(SCALE_CUT_PATH, "test_files/twoPoint_PneE+PeeE.fits")

    BOSS_data_files = [os.path.join(BOSS_PATH,
                                    "CosmoMC_BOSS/data/BOSS.DR12.lowz.3xiwedges_measurements.txt"),
                       os.path.join(BOSS_PATH,
                                    "CosmoMC_BOSS/data/BOSS.DR12.highz.3xiwedges_measurements.txt")]
    BOSS_cov_files = [os.path.join(BOSS_PATH,
                                "CosmoMC_BOSS/data/BOSS.DR12.lowz.3xiwedges_covmat.txt"), 
                      os.path.join(BOSS_PATH,
                                "CosmoMC_BOSS/data/BOSS.DR12.highz.3xiwedges_covmat.txt")]

    KiDS_mock_output_file = os.path.join(mock_dir, "KV450_2x2pt_mock")
    BOSS_mock_output_files = {1 : os.path.join(mock_dir, "BOSS_mock_bin_1.txt"),
                              2 : os.path.join(mock_dir, "BOSS_mock_bin_2.txt")}

    if not create_mocks:
        BOSS_data_files = list(BOSS_mock_output_files.values())
        KiDS_data_cov_file = KiDS_mock_output_file + ".fits"

    source_nz_files = [os.path.join(KCAP_PATH, "data/KV450/nofz/Nz_DIR_z0.1t0.3.asc"),
                        os.path.join(KCAP_PATH, "data/KV450/nofz/Nz_DIR_z0.3t0.5.asc"),
                        os.path.join(KCAP_PATH, "data/KV450/nofz/Nz_DIR_z0.5t0.7.asc"),
                        os.path.join(KCAP_PATH, "data/KV450/nofz/Nz_DIR_z0.7t0.9.asc"),
                        os.path.join(KCAP_PATH, "data/KV450/nofz/Nz_DIR_z0.9t1.2.asc")]

    lens_nz_files = [os.path.join(KCAP_PATH, "data/BOSS/nofz/nOfZ_hist_BOSSA_tomo0.dat"),
                     os.path.join(KCAP_PATH, "data/BOSS/nofz/nOfZ_hist_BOSSA_tomo1.dat")]

    dz_covariance_file = os.path.join(KCAP_PATH, "data/KV450/nofz/DIR_cov.asc")

    source_nz_sample = "KV450_5bin"
    lens_nz_sample = "BOSS"

    BOSS_like_name = "xi_wedges_like"
    if create_mocks:
        bands_range = [0, 160]
        points_range = [0, 32]
    else:
        bands_range = [20, 160]
        points_range = [4, 32]
    num_ell = 3
    z_eff = [0.38, 0.61]

    s = (2.5 + 5*np.arange(60))[points_range[0]:points_range[1]]


    xi_theta_min = 0.5 # arcmin
    xi_theta_max = 300.0 # arcmin
    xi_n_theta = 9

    bandpower_ell_min = 100.0
    bandpower_ell_max = 1500.0
    bandpower_n_bin = 8
    bandpower_theta_min = 0.25
    bandpower_theta_max = 397.0
    bandpower_apodise = 0
    bandpower_Delta_x = 0.5

    used_stats = ["PeeE", "PneE"]

    halofit_version = args.halofit_version

    sampler = args.sampler

    p = K1000Pipeline(KiDS_data_cov_file=KiDS_data_cov_file,
                        BOSS_data_files=BOSS_data_files,
                        BOSS_cov_files=BOSS_cov_files,
                        KiDS_mock_output_file=KiDS_mock_output_file,
                        BOSS_mock_output_files=BOSS_mock_output_files,
                        source_nz_files=source_nz_files,
                        lens_nz_files=lens_nz_files,
                        dz_covariance_file=dz_covariance_file,
                        source_nz_sample=source_nz_sample,
                        lens_nz_sample=lens_nz_sample,
                        BOSS_like_name=BOSS_like_name,
                        bands_range=bands_range,
                        points_range=points_range,
                        num_ell=num_ell,
                        z_eff=z_eff,
                        xi_theta_min=xi_theta_min,
                        xi_theta_max=xi_theta_max,
                        xi_n_theta=xi_n_theta,
                        bandpower_ell_min=bandpower_ell_min,
                        bandpower_ell_max=bandpower_ell_max,
                        bandpower_n_bin=bandpower_n_bin,
                        bandpower_theta_min=bandpower_theta_min,
                        bandpower_theta_max=bandpower_theta_max,
                        bandpower_apodise=bandpower_apodise,
                        bandpower_Delta_x=bandpower_Delta_x,
                        used_stats=used_stats,
                        halofit_version=halofit_version,
                        create_mocks=create_mocks,
                        noisy_mocks=noisy_mocks)
    
    if create_mocks:
        os.makedirs(mock_dir)
        p.choose_pipeline(pipeline_name)
        block = p.run_pipeline()
        if "BOSS_like" in p.config:
            mock_BOSS = p.create_mock_BOSS(block, BOSS_like_name, num_ell, noisy_mocks)
            for b, data in mock_BOSS.items():
                create_BOSS_data_file(data, BOSS_mock_output_files[b])
            for f in BOSS_cov_files:
                shutil.copy(f, mock_dir)


    else:
        p.choose_pipeline(pipeline_name)

        config_dir = os.path.join(run_dir, run_name, "config")
        output_dir = os.path.join(run_dir, run_name, "output")
        print(f"Creating run directory for {run_name}: {config_dir}")
        os.makedirs(config_dir)
        os.makedirs(output_dir)

        config_file = os.path.join(config_dir, f"{run_name}.ini")
        values_file = os.path.join(config_dir, f"{run_name}_values.ini")
        priors_file = os.path.join(config_dir, f"{run_name}_priors.ini")

        derived_parameters=["cosmological_parameters/S_8",
                            "cosmological_parameters/sigma_8",
                            "cosmological_parameters/omega_m",
                            "cosmological_parameters/cosmomc_theta",]
        if "correlated_dz_priors" in p.config:
            # nofz/bias_* values
            derived_parameters += p.config["correlated_dz_priors"]["output_parameters"].split(" ")

        p.set_sampling_config(#likelihoods=["BOSS_like"],
                            derived_parameters=derived_parameters,
                            parameter_file=values_file,
                            prior_file=priors_file,
                            verbose=True,
                            debug=False,
                            output_dir=output_dir,
                            run_name=run_name,
                            sampler_name=sampler,
                            max_iterations="1000000",)
        if sampler == "multinest":
            os.makedirs(os.path.join(output_dir, "multinest"))
    
        ini, values_ini, priors_ini = p.build_ini()

        with open(config_file, "w") as f:
            ini.write(f)
        with open(values_file, "w") as f:
            values_ini.write(f)
        with open(priors_file, "w") as f:
            priors_ini.write(f)
