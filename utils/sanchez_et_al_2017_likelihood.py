from cosmosis.datablock import option_section, names

import numpy as np
import scipy.interpolate

DV_FAP_fsig_data = {0.38 : np.array([9.89, 0.413, 0.468]),
                    0.61 : np.array([14.51, 0.742, 0.439])}

DV_FAP_fsig_cov = { 0.38: np.array([[2.309e-2, -2.201e-4, 8.840e-4], 
                                    [-2.201e-4, 2.005e-4, 4.827e-4],
                                    [8.840e-4, 4.827e-4, 2.763e-3]]),
                    0.61 : np.array([[ 4.253e-2, -9.324e-4, 1.313e-3], 
                                     [-9.324e-4,  5.626e-4, 4.609e-4],
                                     [ 1.313e-3,  4.609e-4, 1.516e-3]])} 

def setup(options):
    mode = options.get_string(option_section, "data_file_name", default="DV_FAP_fsigma8")
    like_name = options.get_string(option_section, "like_name")

    if mode == "DV_FAP_fsigma8":
        inv_cov = {z_ : np.linalg.inv(c) for z_, c in DV_FAP_fsig_cov.items()}
        return DV_FAP_fsig_data, inv_cov, like_name, mode
    elif mode == "DM_H_fsigma8":
        file_name = options.get_string(option_section, "data_file_name", default="FS.txt")
        cov_file_names = [s for s in options.get_string(option_section, "cov_file_names", default="cov1.txt cov.txt").split(" ")]
        z = [float(z) for z in options.get_string(option_section, "z", default="0.38 0.61").split(" ")]

        data  = {}
        with open(file_name) as f:
            for l in f.readlines():
                if l.lstrip()[0] == "#":
                    continue
                z, kind, val = l.split(" ")
                z, val = float(z), float(val)
                if z not in data: data[z] = {}
                data[z].update({kind : val})

        cov = [np.loadtxt(c) for c in cov_file_names]

        if len(cov) != len(z):
            raise ValueError("Mismatch in number of covariances and z values.")

        data = {k : data[k] for k in z}
        inv_cov = {z_ : np.linalg.inv(c) for z_, c in zip(z, cov)}
        return data, inv_cov, like_name, mode
    else:
        raise ValueError(f"Mode {mode} unkown.")

def execute(block, config):
    data_dict, inv_cov_dict, like_name, mode = config

    chi2 = 0
    if mode == "DV_FAP_fsigma8":
        for (z, d), inv_cov in zip(data_dict.items(), inv_cov_dict.values()):
            def interpolate(x_arr, y_arr, x):
                return scipy.interpolate.InterpolatedUnivariateSpline(x_arr, y_arr)(x)

            rs_DV = interpolate(block[names.distances, "z"][1:], block[names.distances, "rs_DV"][1:], z)
            F_AP = interpolate(block[names.distances, "z"][1:], block[names.distances, "F_AP"][1:], z)
            fsigma8 = interpolate(block[names.growth_parameters, "z"], block[names.growth_parameters, "fsigma_8"], z)
            mu = np.array([1/rs_DV, F_AP, fsigma8])

            # print("z:", z, "d:", d, "mu:", mu, "inv_cov:", inv_cov)
            r = mu - d

            chi2 += r @ inv_cov @ r
    else:
        raise ValueError(f"Mode {mode} not supported.")

    ln_like = -0.5*chi2
    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = ln_like

    return 0

def cleanup(config):
    pass