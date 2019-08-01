import numpy as np
from cosmosis.datablock import option_section, names

import p_gm_fit

def setup(options):
    b2_coeffs = np.loadtxt(options[option_section, "b2_coefficient_file"])
    g2_coeffs = np.loadtxt(options[option_section, "g2_coefficient_file"])
    g3_coeffs = np.loadtxt(options[option_section, "g3_coefficient_file"])

    coeffs = {"b2" : b2_coeffs,
              "g2" : g2_coeffs,
              "g3" : g3_coeffs}

    pofk_lin_section = options.get_string("linear_matter_power_section", names.matter_power_lin)
    pofk_nonlin_section = options.get_string("nonlinear_matter_power_section", names.matter_power_nl)
    output_section = options.get_string("output_section", names.matter_galaxy_power)
    
    config = coeffs, pofk_lin_section, pofk_nonlin_section, output_section
    return config

def execute(block, config):
    coeffs, pofk_lin_section, pofk_nonlin_section, output_section = config

    k = block[pofk_lin_section, "k_h"]
    pofk_lin = block[pofk_lin_section, "p_k"]
    pofk_nonlin = block[pofk_nonlin_section, "p_k"]

    omch2 = block[names.cosmological_parameters, "omch2"]
    h = block[names.cosmological_parameters, "h0"]
    n_s = block[names.cosmological_parameters, "n_s"]

    b1 = block["bias_parameters", "b1"]
    b2 = block["bias_parameters", "b2"]
    g2 = block["bias_parameters", "g2"]
    g3 = block["bias_parameters", "g3"]

    P_gm = p_gm_fit.P_gm_approx(k, pofk_lin, pofk_nonlin, omch2, h, n_s, b1, b2, g2, g3, coeffs)
    block[output_section, "p_k"] = P_gm
    
    return 0

def cleanup(config):
    pass
