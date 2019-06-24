from time import sleep

import numpy as np

from cosmosis.datablock import option_section, names
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    n_dim_slow = options.get_int(option_section, "n_dim_slow", default=1)
    parameter_section = options.get_string(option_section, "parameter_section", default="params")
    output_section = options.get_string(option_section, "output_section", default="mu")
    return n_dim_slow, parameter_section, output_section

def execute(block, config):
    n_dim_slow, parameter_section, output_section = config
    theta = [block[parameter_section, f"param_{i}"] for i in range(n_dim_slow)]
    
    sleep(0.001)
    block[output_section, "mu_slow"] = theta
    return 0

def clean(config):
    pass
