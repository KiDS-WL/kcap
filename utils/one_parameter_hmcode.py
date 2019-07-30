import numpy as np
from cosmosis.datablock import option_section, names

def setup(options):
    a_0 = options.get_double(option_section, "a_0", default=0.98)
    a_1 = options.get_double(option_section, "a_1", default=-0.12)
    config = [a_0, a_1]
    print(f"Computing HMCode eta0 as {a_0} + {a_1}*A")
    return config

def execute(block, config):
    a_0, a_1 = config
    A = block[names.halo_model_parameters, "A"]
    eta = a_0 + a_1*A
    block[names.halo_model_parameters, "eta0"] = eta
    return 0

def cleanup(config):
    pass
