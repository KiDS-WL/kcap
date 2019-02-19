import numpy as np

from cosmosis.datablock import option_section

def setup(options):
    alpha = options.get_double(option_section, "alpha", default=0.5)
    return alpha

def execute(block, config):
    alpha = config
    S_star = block["cosmological_parameters", "S_star"]
    Omega_m = block["cosmological_parameters", "omega_m"]
    A_s = S_star*2e-9/(Omega_m/0.3)**alpha
    block["cosmological_parameters", "A_s"] = A_s
    return 0

def clean(config):
    pass
