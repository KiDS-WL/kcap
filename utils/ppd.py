import numpy as np
from cosmosis.datablock import option_section, names

def setup(options):
    like_name = options[option_section, "like_name"]
    return like_name

def execute(block, config):
    like_name = config
    x = block[names.data_vector, like_name + "_theory"]
    s = block[names.data_vector, like_name + "_simulation"]
    inv_cov = block[names.data_vector, like_name + "_inverse_covariance"]
    D = float(np.einsum("i,ij,j", x, inv_cov, s).squeeze())
    block["PPD", "D"] = D
    return 0

def cleanup(config):
    pass
