import numpy as np
import scipy.interpolate

from cosmosis.datablock import option_section
from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    config = {}
    config["Pk_gm_section"] = options.get_string(option_section, "Pk_gm_section")
    config["Pk_mm_lin_section"] = options.get_string(option_section, "Pk_mm_lin_section")
    config["output_section"] = options.get_string(option_section, "output_section")

    config["mode"] = options.get_string(option_section, "mode").lower()
    config["z_interpolation"] = options.get_string(option_section, "z_interpolation").lower()
    config["k_interpolation"] = options.get_string(option_section, "k_interpolation").lower()

    if config["mode"] not in ("pk_gm", "bias", "linear_bias"):
        raise ValueError(f"Mode {config['mode']} not supported.")
    if config["k_interpolation"] not in ("none", "const", "lin"):
        raise ValueError(f"Mode {config['k_interpolation']} not supported.")
    if config["z_interpolation"] not in ("none", "const", "lin"):
        raise ValueError(f"Mode {config['z_interpolation']} not supported.")

    if config["mode"] in ("bias", "linear_bias"):
        config["Pk_mm_section"] = options.get_string(option_section, "Pk_mm_section")
        config["Pk_mm_nonlin_section"] = options.get_string(option_section, "Pk_mm_nonlin_section")
    

    n_bin = options[option_section, "n_bin"]
    if config["z_interpolation"] == "const" and n_bin != 2:
        raise ValueError("z_interpolation mode const only works for 2 input redshifts.")

    z_range = {}
    for i in range(n_bin):
        z_range[i+1] = options[option_section, f"z_range_bin_{i+1}"]

    config["z_ranges"] = z_range

    return config

def execute(block, config):
    z_gm = block[config["Pk_gm_section"], "z"]
    k_h_gm = block[config["Pk_gm_section"], "k_h"]
    p_k_gm = block[config["Pk_gm_section"], "p_k"]

    p_k_gm[p_k_gm<0] = 0

    if config["mode"] == "pk_gm":
        z_intp = block[config["Pk_mm_lin_section"], "z"]

        block[config["output_section"], "z"] = z_intp
        block[config["output_section"], "k_h"] = k_h_gm

        if config["z_interpolation"] == "none":
            p_k_gm_extrapolate = np.zeros((len(z_intp), len(k_h_gm)))
            for i, (z_min, z_max) in enumerate(config["z_ranges"].values()):
                p_k_gm_extrapolate[(z_min <= z_intp) & (z_intp < z_max)] = p_k_gm[i]
            block[config["output_section"], "p_k"] = p_k_gm_extrapolate
        elif config["z_interpolation"] == "const":
            p_k_gm_extrapolate = np.zeros((len(z_intp), len(k_h_gm)))
            p_k_gm_extrapolate[z_intp < config["z_ranges"][1][1]] = p_k_gm[0]
            p_k_gm_extrapolate[config["z_ranges"][2][0] <= z_intp] = p_k_gm[1]
            block[config["output_section"], "p_k"] = p_k_gm_extrapolate
        elif config["z_interpolation"] == "lin":
            intp = scipy.interpolate.RectBivariateSpline(z_gm, k_h_gm, p_k_gm, kx=1, ky=1)
            block[config["output_section"], "p_k"] = intp(z_intp, k_h_gm, grid=True)

    elif config["mode"] == "bias":
        z_intp = block[config["Pk_mm_nonlin_section"], "z"]
        k_h_mm_nonlin = block[config["Pk_mm_nonlin_section"], "k_h"]
        p_k_mm_nonlin = block[config["Pk_mm_nonlin_section"], "p_k"]

        p_k_mm = block[config["Pk_mm_section"], "p_k"]
        p_k_mm[p_k_mm<0] = 0
        bias = p_k_gm/p_k_mm

        bounds = np.argwhere(np.diff(np.isclose(p_k_mm,0))).squeeze()
        # print("Support of bias:", bounds)
        # print("Values at bounds:")
        for i in range(len(z_gm)):
            # print("Lower:", bias[i,bounds[2*i][1]-1], bias[i,bounds[2*i][1]], bias[i,bounds[2*i][1]+1])
            # print("Upper:", bias[i,bounds[2*i+1][1]-1], bias[i,bounds[2*i+1][1]], bias[i,bounds[2*i+1][1]+1])
            if config["k_interpolation"] == "none":
                #Lower 
                bias[i,:bounds[2*i][1]+1] = 0.0
                #Upper
                bias[i,bounds[2*i+1][1]+1:] = 0.0
            elif config["k_interpolation"] == "const":
                #Lower 
                bias[i,:bounds[2*i][1]+1] = bias[i,bounds[2*i][1]+1]
                #Upper
                bias[i,bounds[2*i+1][1]+1:] = bias[i,bounds[2*i+1][1]]

        if config["z_interpolation"] == "none":
            bias_extrapolate = np.zeros((len(z_intp), len(k_h_gm)))
            for i, (z_min, z_max) in enumerate(config["z_ranges"].values()):
                bias_extrapolate[(z_min <= z_intp) & (z_intp < z_max)] = bias[i]
        elif config["z_interpolation"] == "const":
            bias_extrapolate = np.zeros((len(z_intp), len(k_h_gm)))
            bias_extrapolate[z_intp < config["z_ranges"][1][1]] = bias[0]
            bias_extrapolate[config["z_ranges"][2][0] <= z_intp] = bias[1]
        elif config["z_interpolation"] == "lin":
            intp = scipy.interpolate.RectBivariateSpline(z_gm, k_h_gm, bias, kx=1, ky=1)
            bias_extrapolate = intp(z_intp, k_h_gm, grid=True)

        block[config["output_section"], "z"] = z_intp
        block[config["output_section"], "k_h"] = k_h_mm_nonlin
        block[config["output_section"], "p_k"] = bias_extrapolate*p_k_mm_nonlin

    elif config["mode"] == "linear_bias":
        b = block["bias_parameters", "b1"]
        block[config["output_section"], "z"] = block[config["Pk_mm_nonlin_section"], "z"]
        block[config["output_section"], "k_h"] = block[config["Pk_mm_nonlin_section"], "k_h"]
        block[config["output_section"], "p_k"] = b*block[config["Pk_mm_nonlin_section"], "p_k"]
        
    block[config["output_section"], "_cosmosis_order_p_k"] = "z_cosmosis_order_k_h"

    return 0

def cleanup(config):
    pass
