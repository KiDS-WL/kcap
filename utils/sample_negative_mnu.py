import numpy as np

def setup(options):
    config = []
    return config

def execute(block, config):
    mnu = block["cosmological_parameters", "mnu_proxy"]
    block["cosmological_parameters", "mnu"] = np.abs(mnu)
    return 0

def cleanup(config):
    pass
