import numpy as np

from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    constant_c_offset = options.get_logical_default(option_section, "constant_c_offset", "T")

    try:
        filename = options.get_string(option_section, "xi_pm_c_file")
        xi_p, xi_m = np.loadtxt(filename, unpack=True, usecols=[3,4])
    except errors.BlockNameNotFound:
        xi_p = None
        xi_m = None

    config = constant_c_offset, xi_p, xi_m
    return config

def execute(block, config):
    constant_c_offset, xi_p_c, xi_m_c = config

    n_bin = block["shear_xi", "nbin_a"]

    if constant_c_offset:
        delta_c = block["shear_c_bias", "delta_c"]
    else:
        delta_c = 0.0

    if xi_p_c is not None:
        A_c = block["shear_c_bias", "A_c"]
    else:
        A_c = 0.0

    for i in range(n_bin):
        for j in range(i):
            block["shear_xi", f"xiplus_{i+1}_{j+1}"] += A_c**2 * xi_p_c + delta_c
            block["shear_xi", f"ximinus_{i+1}_{j+1}"] += A_c**2 * xi_m_c + delta_c

    return 0

def cleanup(config):
    pass
