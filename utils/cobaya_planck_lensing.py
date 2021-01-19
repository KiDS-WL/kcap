import numpy as np

from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

# This uses the cobaya standalone CMB lensing likelihood https://github.com/CobayaSampler/planck_lensing_external
import plancklensing

cosmo = section_names.cosmological_parameters

def setup(options):
    like_name = options.get_string(option_section, "like_name", default="Planck2018_lensing")

    likelihood = plancklensing.PlanckLensingMarged()
    return likelihood, like_name


def execute(block, config):
    likelihood, like_name = config

    cl_PP = block[section_names.cmb_cl, "pp"]
    cl_PP = np.concatenate(([0,0], cl_PP))

    ell = np.arange(len(cl_PP))

    logL = likelihood.log_likelihood({"pp" : ell*(ell+1)*cl_PP}, A_Planck=1.0)

    block[section_names.data_vector, like_name+"_CHI2"] = -2*logL
    block[section_names.likelihoods, like_name + "_LIKE"] = logL 
    
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
