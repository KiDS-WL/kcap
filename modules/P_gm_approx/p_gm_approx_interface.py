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

    #pofk_lin_section = options.get_string("linear_matter_power_section", names.matter_power_lin)
    pofk_lin_section = options.get_string(option_section,"linear_matter_power_section")
    #pofk_nonlin_section = options.get_string("nonlinear_matter_power_section", names.matter_power_nl)
    pofk_nonlin_section = options.get_string(option_section,"nonlinear_matter_power_section")
    #output_section = options.get_string("output_section", names.matter_galaxy_power)
    output_section = options.get_string(option_section,"output_section")
    
    # read flag for Lagrangian bias approx.
    flag_lag_g2 = options.get_bool(option_section, "local_lag_g2", True)

    # read redshift boundary between low-z and high-z bin
    z_sep = options.get_double(option_section, "z_sep")

    config = coeffs, pofk_lin_section, pofk_nonlin_section, output_section, flag_lag_g2, z_sep
    return config

def execute(block, config):
    coeffs, pofk_lin_section, pofk_nonlin_section, output_section, flag_lag_g2, z_sep = config

    z = block[pofk_nonlin_section, "z"]
    k = block[pofk_nonlin_section, "k_h"]
    z_lin = block[pofk_lin_section, "z"]
    pofk_lin = block[pofk_lin_section, "p_k"]
    pofk_nonlin = block[pofk_nonlin_section, "p_k"]
    if (pofk_lin.shape != pofk_nonlin.shape):  # should only happen when using PT non-lin. power spectrum
        pofk_lin = p_gm_fit.pofk_interpolator(pofk_lin,k,z_lin)(k,z)  # put lin. PS on same grid as nonlin. PS

    omch2 = block[names.cosmological_parameters, "omch2"]
    h = block[names.cosmological_parameters, "h0"]
    n_s = block[names.cosmological_parameters, "n_s"]

    b1_lowz = block["bias_parameters", "b1_bin_1"]
    b1_higz = block["bias_parameters", "b1_bin_2"]

    b2_lowz = block["bias_parameters", "b2_bin_1"]
    b2_higz = block["bias_parameters", "b2_bin_2"]

    g3_lowz = block["bias_parameters", "gamma3_bin_1"]
    g3_higz = block["bias_parameters", "gamma3_bin_2"]

    if (flag_lag_g2):
        g2_lowz = -2./7.*(b1_lowz-1.)
        g2_higz = -2./7.*(b1_higz-1.)
    else:
        g2_lowz = block["bias_parameters", "gamma2_bin_1"]
        g2_higz = block["bias_parameters", "gamma2_bin_2"]

    P_gm_lowz = p_gm_fit.P_gm_approx(k, pofk_lin, pofk_nonlin, omch2, h, n_s, b1_lowz, b2_lowz, g2_lowz, g3_lowz, coeffs)
    P_gm_higz = p_gm_fit.P_gm_approx(k, pofk_lin, pofk_nonlin, omch2, h, n_s, b1_higz, b2_higz, g2_higz, g3_higz, coeffs)

    # merge such that returned array contains low-z power spectrum below z_sep and high-z power spectrum above z_sep
    zcond = z<z_sep
    P_gm = np.where(zcond[:,np.newaxis],P_gm_lowz,P_gm_higz)

    block[output_section, "p_k"] = P_gm

    block[output_section, "z"] = z
    block[output_section, "k_h"] = k
    
    return 0

def cleanup(config):
    pass
