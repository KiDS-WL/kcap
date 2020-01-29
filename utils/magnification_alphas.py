import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    alphas = options[option_section, "alpha_binned"]
        
    return alphas

def execute(block, config):
    alphas = config
    block[names.galaxy_luminosity_function, "alpha_binned"] = alphas
    return 0

def clean(config):
    pass
