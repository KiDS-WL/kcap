import numpy as np

def setup(options):
    config = []
    return config

def execute(block, config):
    ln_1e10_A_s = block["cosmological_parameters", "ln_1e10_A_s"]
    A_s = np.exp(ln_1e10_A_s)*1e-10
    block["cosmological_parameters", "A_s"] = A_s
    return 0

def cleanup(config):
    pass
