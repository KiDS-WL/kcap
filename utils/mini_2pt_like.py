import collections

import scipy.interpolate
import numpy as np
pi = np.pi

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    n_bin = options[option_section, "n_bin"]
    try:
        single_filename = options[option_section, "data_filename"]
        data = np.loadtxt(single_filename).T
        theta = data[0]
        n_theta = len(theta)//2
        theta_xi_plus = theta[:n_theta]*pi/180/60
        theta_xi_minus = theta[n_theta:]*pi/180/60
    except errors.BlockNameNotFound:
        single_filename = None

    try:
        cut_xi_plus_idx = options[option_section, "cut_xi_plus"]
    except errors.BlockNameNotFound:
        cut_xi_plus_idx = []
    try:
        cut_xi_minus_idx = options[option_section, "cut_xi_minus"]
    except errors.BlockNameNotFound:
        cut_xi_minus_idx = []

    order_data = options.get_string(option_section, "order_data", "xi_pm-bin-theta").lower()
    order_cov = options.get_string(option_section, "order_cov", "xi_pm-bin-theta").lower()

    cov_filename = options.get_string(option_section, "cov")
    cov = np.loadtxt(cov_filename)
    if order_cov == "xi_pm-bin-theta":
        if cut_xi_minus_idx != []:
            idx_xi_minus = (n_bin*(n_bin+1))//2*n_theta + np.arange((n_bin*(n_bin+1))//2, dtype=int)*n_theta + cut_xi_minus_idx[:,None]
        else:
            idx_xi_minus = np.array([], ndmin=1)
        if cut_xi_plus_idx != []:
            idx_xi_plus = np.arange((n_bin*(n_bin+1))//2, dtype=int)*n_theta + cut_xi_plus_idx[:,None]
        else:
            idx_xi_plus = np.array([], ndmin=1)
        idx = np.concatenate([idx_xi_plus.flatten(), idx_xi_minus.flatten()])
        cov = np.delete(cov, idx, axis=0)
        cov = np.delete(cov, idx, axis=1)
    elif order_cov == "montepython":
        if cut_xi_minus_idx != []:
            idx_xi_minus = np.arange((n_bin*(n_bin+1))//2, dtype=int)*2*n_theta + n_theta + cut_xi_minus_idx[:,None]
        else:
            idx_xi_minus = np.array([], ndmin=1)
        if cut_xi_plus_idx != []:
            idx_xi_plus = np.arange((n_bin*(n_bin+1))//2, dtype=int)*2*n_theta + cut_xi_plus_idx[:,None]
        else:
            idx_xi_plus = np.array([], ndmin=1)
        idx = np.concatenate([idx_xi_plus.flatten(), idx_xi_minus.flatten()])
        cov = np.delete(cov, idx, axis=0)
        cov = np.delete(cov, idx, axis=1)
    else:
        raise ValueError(f"Unsupported covariance order {order_cov}.")
    
    inv_cov = np.linalg.inv(cov)

    data_vectors = {}
    for i in range(n_bin):
        data_vectors[i] = {}
        for j in range(i+1):
            if single_filename is not None:
                if order_data == "xi_pm-bin-theta":
                    bin_idx =  (i*(i+1))//2+j
                elif order_data == "montepython":
                    bin_idx = (j*(2*n_bin - j + 1))//2 + i-j
                else:
                    raise ValueError(f"Unsupported data order {order_data}.")
                xi_plus = data[1+bin_idx,:n_theta]
                xi_minus = data[1+bin_idx,n_theta:]
                data_vectors[i][j] = [xi_plus, xi_minus]
            else:
                filename = options.get_string(option_section, f"bin_{i+1}_{j+1}")
                data = np.loadtxt(filename).T
                theta = data[0]*pi/180/60
                theta_xi_plus = theta
                theta_xi_minus = theta
                data_vectors[i][j] = data[[1,2]]

            data_vectors[i][j][0] = np.delete(data_vectors[i][j][0], cut_xi_plus_idx)
            data_vectors[i][j][1] = np.delete(data_vectors[i][j][1], cut_xi_minus_idx)
                
    theta_xi_plus = np.delete(theta_xi_plus, cut_xi_plus_idx)
    theta_xi_minus = np.delete(theta_xi_minus, cut_xi_minus_idx)

    like_name = options.get_string(option_section, "like_name")
    keep_theory_vector = options.get_bool(option_section, "keep_theory_vector", False)
    return inv_cov, data_vectors, theta_xi_plus, theta_xi_minus, order_cov, like_name, keep_theory_vector

def execute(block, config):
    inv_cov, data_vectors, theta_xi_plus, theta_xi_minus, order_cov, like_name, keep_theory_vector = config
    
    n_bin = block["shear_xi", "nbin_a"]

    if order_cov == "xi_pm-bin-theta":
        theory_xi_plus_vector = []
        theory_xi_minus_vector = []
        data_xi_plus_vector = []
        data_xi_minus_vector = []
        
        for i in range(n_bin):
            for j in range(i+1):
                data_xi_plus, data_xi_minus = data_vectors[i][j]

                theory_theta = block["shear_xi", "theta"]
                theory_xi_plus = block["shear_xi", f"xiplus_{i+1}_{j+1}"]
                theory_xi_minus = block["shear_xi", f"ximinus_{i+1}_{j+1}"]

                intp_xi_plus = scipy.interpolate.InterpolatedUnivariateSpline(np.log(theory_theta), theory_xi_plus)
                theory_xi_plus = intp_xi_plus(np.log(theta_xi_plus))
                intp_xi_minus = scipy.interpolate.InterpolatedUnivariateSpline(np.log(theory_theta), theory_xi_minus)
                theory_xi_minus = intp_xi_minus(np.log(theta_xi_minus))

                data_xi_plus_vector.append(data_xi_plus)
                data_xi_minus_vector.append(data_xi_minus)
                theory_xi_plus_vector.append(theory_xi_plus)
                theory_xi_minus_vector.append(theory_xi_minus)
                    
        x = np.concatenate(data_xi_plus_vector + data_xi_minus_vector)
        mu = np.concatenate(theory_xi_plus_vector + theory_xi_minus_vector)
    elif order_cov == "montepython":
        data_xi_vector = []
        theory_xi_vector = []
        for i in range(n_bin):
            for j in range(i,n_bin):
                data_xi_plus, data_xi_minus = data_vectors[j][i]

                theory_theta = block["shear_xi", "theta"]
                theory_xi_plus = block["shear_xi", f"xiplus_{j+1}_{i+1}"]
                theory_xi_minus = block["shear_xi", f"ximinus_{j+1}_{i+1}"]

                intp_xi_plus = scipy.interpolate.InterpolatedUnivariateSpline(np.log(theory_theta), theory_xi_plus)
                theory_xi_plus = intp_xi_plus(np.log(theta_xi_plus))
                intp_xi_minus = scipy.interpolate.InterpolatedUnivariateSpline(np.log(theory_theta), theory_xi_minus)
                theory_xi_minus = intp_xi_minus(np.log(theta_xi_minus))

                data_xi_vector.append(np.concatenate([data_xi_plus, data_xi_minus]))
                theory_xi_vector.append(np.concatenate([theory_xi_plus, theory_xi_minus]))
        
        x = np.concatenate(data_xi_vector)
        mu = np.concatenate(theory_xi_vector)

    d = x - mu

    chi2 = np.einsum('i,ij,j', d, inv_cov, d)
    chi2 = float(chi2)
    like = -0.5*chi2

    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = like

    if keep_theory_vector:
        block[names.data_vector, like_name + "_theory"] = mu
        block[names.data_vector, like_name + "_data"] = x
        block[names.data_vector, like_name + "_theta_xi_plus"] = theta_xi_plus
        block[names.data_vector, like_name + "_theta_xi_minus"] = theta_xi_minus

    return 0

def cleanup(config):
    pass
