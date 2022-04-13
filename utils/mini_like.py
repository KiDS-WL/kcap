import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    like_name = options.get_string(option_section, "like_name")
    input_section_name = options.get_string(option_section, "input_section_name", default="likelihood")
        
    return like_name, input_section_name

def execute(block, config):
    like_name, input_section_name = config

    d = block[input_section_name, "data"]
    mu = block[input_section_name, "theory"]

    inv_cov = block[input_section_name, "inv_covariance"]
    r = d - mu

    chi2 = float(r @ inv_cov @ r)
    ln_like = -0.5*chi2

    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.data_vector, like_name+"_THEORY"] = mu
    block[names.data_vector, like_name+"_INVERSE_COVARIANCE"] = inv_cov
    block[names.likelihoods, like_name+"_LIKE"] = ln_like

    return 0

def clean(config):
    pass
