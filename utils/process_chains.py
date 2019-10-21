import re

import numpy as np

parameter_dictionary = {
        "omega_m" :                  {"cosmosis" :    "cosmological_parameters--omega_m",
                                      "montepython" : "Omega_m",
                                      "cosmomc" :     "omegam",
                                      "latex" :       "\\Omega_m"},
        "omega_c h^2" :              {"cosmosis" :    "cosmological_parameters--omch2",
                                      "montepython" : "omega_cdm",
                                      "cosmomc" :     "omegach2",
                                      "latex" :       "h^2\\Omega_c"},
        "omega_b h^2" :              {"cosmosis" :    "cosmological_parameters--ombh2", 
                                      "montepython" : "omega_b",
                                      "cosmomc" :     "omegabh2",
                                      "latex" :       "h^2\\Omega_b"},
        "h" :                        {"cosmosis" :    "cosmological_parameters--h0", 
                                      "montepython" : "h",
                                      "cosmomc" :     "h",
                                      "latex" :       "h"},
        "w" :                        {"cosmosis" :    "cosmological_parameters--w", 
                                      "montepython" : "w",
                                      "cosmomc" :     "w",
                                      "latex" :       "w"},
        "w_background" :             {"cosmosis" :    "cosmological_parameters--w_background", 
                                      "montepython" : "w",
                                      "cosmomc" :     "w_b",
                                      "latex" :       "w_b"},
        "n_s"   :                    {"cosmosis" :    "cosmological_parameters--n_s",
                                      "montepython" : "n_s",
                                      "cosmomc" :     "ns",
                                      "latex" :       "n_s"},
        "A_s"   :                    {"cosmosis" :    "cosmological_parameters--a_s",
                                      "montepython" : "A_s",
                                      "cosmomc" :     "as",
                                      "latex" :       "A_s"},
        "log A_s"   :                {"cosmosis" :    "cosmological_parameters--ln_1e10_a_s",
                                      "montepython" : "ln10^{10}A_s",
                                      "cosmomc" :     "logA",
                                      "latex" :       "{\\rm ln}(10^{10} A_s)"},
        "sigma_8" :                  {"cosmosis" :    "cosmological_parameters--sigma_8",
                                      "montepython" : "sigma8",
                                      "cosmomc" :     "sigma8",
                                      "latex" :       "\\sigma_8"},
        "S_8" :                      {"cosmosis" :    "cosmological_parameters--s8",
                                      "montepython" : "S8",
                                      "cosmomc" :     "s8",
                                      "latex" :       "S_8"},
        "cosmomc_theta" :            {"cosmosis" :    "cosmological_parameters--cosmomc_theta",
                                      "cosmomc" :     "theta",
                                      "latex" :       "100\\theta_{\\rm MC}"},
        "tau" :                      {"cosmosis" :    "cosmological_parameters--tau",
                                      "cosmomc" :     "tau",
                                      "latex" :       "\\tau"},
        "b_1 lowz" :                 {"cosmosis" :    "bias_parameters--b1_bin_1",
                                      "cosmomc" :     "b1l",
                                      "latex" :       "b_1^{\\rm lowz}"},
        "b_2 lowz" :                 {"cosmosis" :    "bias_parameters--b2_bin_1",
                                      "cosmomc" :     "b2l",
                                      "latex" :       "b_2^{\\rm lowz}"},
        "gamma_3 lowz" :             {"cosmosis" :    "bias_parameters--gamma3_bin_1",
                                      "cosmomc" :     "gamma3l",
                                      "latex" :       "\\gamma_3^{\\rm lowz}"},
        "a_vir lowz" :               {"cosmosis" :    "bias_parameters--a_vir_bin_1",
                                      "cosmomc" :     "a_virl",
                                      "latex" :       "a_{vir}^{\\rm lowz}"},
        "b_1 highz" :                {"cosmosis" :    "bias_parameters--b1_bin_2",
                                      "cosmomc" :     "b1h",
                                      "latex" :       "b_1^{\\rm highz}"},
        "b_2 highz" :                {"cosmosis" :    "bias_parameters--b2_bin_2",
                                      "cosmomc" :     "b2h",
                                      "latex" :       "b_2^{\\rm highz}"},
        "gamma_3 highz" :            {"cosmosis" :    "bias_parameters--gamma3_bin_2",
                                      "cosmomc" :     "gamma3h",
                                      "latex" :       "\\gamma_3^{\\rm highz}"},
        "a_vir highz" :              {"cosmosis" :    "bias_parameters--a_vir_bin_2",
                                      "cosmomc" :     "a_virh",
                                      "latex" :       "a_{\rm vir}^{\\rm highz}"},
        "A_IA" :                     {"cosmosis" :    "intrinsic_alignment_parameters--a",
                                      "montepython" : "A_IA",
                                      "cosmomc" :     "a_ia",
                                      "latex" :       "A_{\\rm IA}"},
        "A_baryon" :                 {"cosmosis" :    "halo_model_parameters--a",
                                      "montepython" : "c_min",
                                      "cosmomc" :     "a_baryon",
                                      "latex" :       "A_{\\rm baryon}"},
        "eta_baryon" :               {"cosmosis" :    "halo_model_parameters--eta0",
                                      "montepython" : "eta_0",
                                      "cosmomc" :     "eta_baryon",
                                      "latex" :       "\\eta_{\\rm baryon}"},
        "delta_c" :                  {"cosmosis" :    "shear_c_bias--delta_c",
                                      "montepython" : "dc",
                                      "cosmomc" :     "delta_c",
                                      "latex" :       "\\delta c"},
        "A_c" :                      {"cosmosis" :    "shear_c_bias--a_c",
                                      "montepython" : "Ac",
                                      "cosmomc" :     "a_c",
                                      "latex" :       "A_{\\rm c}"},
        "delta z_1" :                {"cosmosis" :    "nofz_shifts--bias_1",
                                      "montepython" : "D_z1",
                                      "cosmomc" :     "delta_z1",
                                      "latex" :       "\\delta z_1"},
        "delta z_2" :                {"cosmosis" :    "nofz_shifts--bias_2",
                                      "montepython" : "D_z2",
                                      "cosmomc" :     "delta_z2",
                                      "latex" :       "\\delta z_2"},
        "delta z_3" :                {"cosmosis" :    "nofz_shifts--bias_3",
                                      "montepython" : "D_z3",
                                      "cosmomc" :     "delta_z3",
                                      "latex" :       "\\delta z_3"},
        "delta z_4" :                {"cosmosis" :    "nofz_shifts--bias_4",
                                      "montepython" : "D_z4",
                                      "cosmomc" :     "delta_z4",
                                      "latex" :       "\\delta z_4"},
        "delta z_5" :                {"cosmosis" :    "nofz_shifts--bias_5",
                                      "montepython" : "D_z5",
                                      "cosmomc" :     "delta_z5",
                                      "latex" :       "\\delta z_5"},
        "A_Planck" :                 {"cosmosis" :    "planck--a_planck",
                                      "cosmomc" :     "calPlanck",
                                      "latex" :       "A_{\\rm Planck}"},
        "F_AP lowz" :                {"cosmosis" :    "lss_parameters--f_ap_bin_1",
                                      "cosmomc" :     "FAP1",
                                      "latex" :       "F_{\\rm AP lowz}"},
        "rs_DV lowz" :               {"cosmosis" :    "lss_parameters--rs_dv_bin_1",
                                      "cosmomc" :     "rsDv1",
                                      "latex" :       "r_{\\rm d}/D_V {\\rm lowz}"},
        "fsigma_8 lowz" :            {"cosmosis" :    "lss_parameters--fsigma_8_bin_1",
                                      "cosmomc" :     "fsigma8z1",
                                      "latex" :       "f\\sigma_8 {\\rm highz}"},
        "F_AP highz" :               {"cosmosis" :    "lss_parameters--f_ap_bin_2",
                                      "cosmomc" :     "FAP2",
                                      "latex" :       "F_{\\rm AP highz}"},
        "rs_DV highz" :              {"cosmosis" :    "lss_parameters--rs_dv_bin_2",
                                      "cosmomc" :     "rsDv2",
                                      "latex" :       "r_{\\rm d}/D_V {\\rm lowz}"},
        "fsigma_8 highz" :           {"cosmosis" :    "lss_parameters--fsigma_8_bin_2",
                                      "cosmomc" :     "fsigma8z2",
                                      "latex" :       "f\\sigma_8 {\\rm highz}"},
           
                              }

def parameter_mapping(chain_file,
                      parameters,
                      statistics=[], chain_format="cosmosis", 
                      parameter_name_map=[]):
    chain_format = chain_format.lower()
    parameter_names = parameter_dictionary
    
    with open(chain_file, "r") as f:
        params = f.readline()[1:]
    
    if chain_format == "cosmosis":
        params = [p.strip().lower() for p in params.split("\t")]
    elif chain_format == "montepython":
        params = [p.strip() for p in params.split(",")]
    else:
        raise ValueError(f"Chain format {chain_format} not supported.")

    param_idx = []
    output_parameter_names = {mapping : [] for mapping in parameter_name_map}
    for p in parameters:
        if p in parameter_names and parameter_names[p][chain_format] in params:
            param_idx.append(params.index(parameter_names[p][chain_format]))
            if parameter_name_map is not None:
                for mapping in parameter_name_map:
                    output_parameter_names[mapping].append(parameter_names[p][mapping])
        else:
            raise ValueError(f"Parameter {p} not in chain.")
            
    stats_idx = []
    for s in statistics:
        s = s.lower()
        if s in params:
            stats_idx.append(params.index(s))
        else:
            raise ValueError(f"Statistic {s} not in chain.")
    
    if parameter_name_map is not []:
        return (param_idx, stats_idx) + tuple(output_parameter_names.values())
    else:
        return param_idx, stats_idx

def load_nested_sampling_file(filename):
    info = {}
    with open(filename, "r") as f:
        for s in f.readlines()[-3:]:
                m = re.match("^#([a-z_]+)=([0-9\.\-]+)", s)
                if m is None:
                    raise ValueError(f"Can't read property from line {s}")
                key, value = m.groups()
                info[key] = value
            
    n_sample = int(info["nsample"])
    chain = np.loadtxt(filename)
    print(f"Using {n_sample} samples out of {chain.shape[0]} in the chain.")
    return chain[-n_sample:], float(info["log_z"]), float(info["log_z_error"])