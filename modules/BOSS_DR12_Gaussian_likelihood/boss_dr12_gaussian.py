import os

from cosmosis.datablock import option_section, names

import numpy as np
import scipy.interpolate

HERE = os.path.split(os.path.abspath(__file__))[0]

DV_FAP_fsig_data = {0.38 : np.array([9.89, 0.413, 0.468]),
                    0.61 : np.array([14.51, 0.742, 0.439])}

DV_FAP_fsig_cov = { 0.38: np.array([[2.309e-2, -2.201e-4, 8.840e-4], 
                                    [-2.201e-4, 2.005e-4, 4.827e-4],
                                    [8.840e-4, 4.827e-4, 2.763e-3]]),
                    0.61 : np.array([[ 4.253e-2, -9.324e-4, 1.313e-3], 
                                     [-9.324e-4,  5.626e-4, 4.609e-4],
                                     [ 1.313e-3,  4.609e-4, 1.516e-3]])} 

def interpolate(x_arr, y_arr, x):
    return scipy.interpolate.InterpolatedUnivariateSpline(x_arr, y_arr)(x)

def setup(options):
    mode = options.get_string(option_section, "mode", default="DV_FAP_fsigma8")
    like_name = options.get_string(option_section, "like_name")

    if mode == "DV_FAP_fsigma8":
        default_data_file = os.path.join(HERE, "data/COMBINEDDR12_final_consensus_dV_FAP/final_consensus_results_dV_FAP_fsig.txt")
        default_cov_file = os.path.join(HERE, "data/COMBINEDDR12_final_consensus_dV_FAP/final_consensus_covtot_dV_FAP_fsig.txt")
    else:
        raise ValueError(f"Mode {mode} not supported.")
    
    data_file = options.get_string(option_section, "data_file", default=default_data_file)
    cov_file = options.get_string(option_section, "cov_file", default=default_cov_file)

    redshifts = [float(z) for z in options.get_string(option_section, "z", default="0.38 0.51 0.61").split(" ")]

    data  = {}
    with open(data_file) as f:
        for l in f.readlines():
            if l.lstrip()[0] == "#":
                continue
            z, kind, val = l.split()
            z, val = float(z), float(val)
            if z not in data: data[z] = {}
            data[z].update({kind : val})

    cov = [np.loadtxt(c) for c in cov_file.split(" ")]

    n_param = 3
    n_z = len(data.keys())

    # Assemble covariances into a single covariance matrix
    if len(cov) == 1 and cov[0].shape == (n_param*n_z, n_param*n_z):
        cov_total = cov[0]
    elif len(cov) > 1 and cov[0].shape == (n_param, n_param):
        cov_total = np.zeros((n_param*n_z, n_param*n_z))
        for i in range(n_z):
            cov_total[i*n_z:(i+1)*n_z,i*n_z:(i+1)*n_z] = cov[i]
    else:
        raise ValueError("Mismatch in number of covariances and z values.")

    data = np.concatenate([list(data[z].values()) for z in redshifts])
    print("data:", data)
    inv_cov = np.linalg.inv(cov_total)

    return data, inv_cov, redshifts, like_name, mode


def execute(block, config):
    data, inv_cov, redshifts, like_name, mode = config

    n_param = len(data) // len(redshifts)

    if mode == "DV_FAP_fsigma8":
        mu = np.zeros_like(data)

        for i, z in enumerate(redshifts):
            rs_DV = interpolate(block[names.distances, "z"][1:], block[names.distances, "rs_DV"][1:], z)
            F_AP = interpolate(block[names.distances, "z"][1:], block[names.distances, "F_AP"][1:], z)
            fsigma8 = interpolate(block[names.growth_parameters, "z"], block[names.growth_parameters, "fsigma_8"], z)
            
            mu[i*n_param:(i+1)*n_param] = np.array([1/rs_DV, F_AP, fsigma8])
    else:
        raise ValueError(f"Mode {mode} not supported.")

    # print("z:", z, "d:", d, "mu:", mu, "inv_cov:", inv_cov)
    r = data - mu
    chi2 = r @ inv_cov @ r

    ln_like = -0.5*chi2
    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = ln_like

    return 0

def cleanup(config):
    pass