import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    try:
        alphas = options[option_section, "alpha_binned"]
    except errors.BlockNameNotFound:
        # Use hard coded values:
        alphas = np.array([1.0, 2.0])
        
    return alphas

def execute(block, config):
    alphas = config
    block[names.galaxy_luminosity_function, "alpha_binned"] = alphas
    return 0

def clean(config):
    pass
