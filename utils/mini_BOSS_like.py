import collections

import scipy.interpolate
import numpy as np
pi = np.pi

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors


def setup(options):
    data_file = options[option_section, "data_vector_file"]
    data_vector = np.loadtxt(data_file)

    cov_file = options[option_section, "covariance_file"]
    cov = np.loadtxt(cov_file)

    try:
        points_range = options[option_section, "points_range"]
    except errors.BlockNameNotFound:
        points_range = [4, 32]
    
    n_wedge = (data_vector.shape[1]-1)//2
    n_points_full = data_vector.shape[0]
    data_vector = np.concatenate([data_vector[:,2*i+1] for i in range(n_wedge)])
    
    cut_points_idx = list(range(points_range[0])) + list(range(points_range[1], n_points_full))

    if len(cut_points_idx) > 0:
        cut_points_idx = np.array(cut_points_idx, dtype=int)
        cut_points_idx = np.concatenate([cut_points_idx + n_points_full*i for i in range(n_wedge)])

        data_vector = np.delete(data_vector, cut_points_idx) 
        cov = np.delete(cov, cut_points_idx, axis=0)
        cov = np.delete(cov, cut_points_idx, axis=1)
    
    inv_cov = np.linalg.inv(cov)
    
    like_name = options.get_string(option_section, "like_name")
    keep_theory_vector = options.get_bool(option_section, "keep_theory_vector", False)
    return inv_cov, data_vector, cut_points_idx, like_name, keep_theory_vector

def execute(block, config):
    inv_cov, data_vector, cut_points_idx, like_name, keep_theory_vector = config
    
    n_wedge = block["xi_wedges", "n_wedge"]
    mu = block["xi_wedges", "vtheo_convolved"]
    mu = np.concatenate([m for m in mu])
    x = data_vector
    d = x - mu

    chi2 = np.einsum('i,ij,j', d, inv_cov, d)
    chi2 = float(chi2)
    like = -0.5*chi2

    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = like

    if keep_theory_vector:
        block[names.data_vector, like_name + "_theory"] = mu
        block[names.data_vector, like_name + "_data"] = x
        block[names.data_vector, like_name + "_inv_cov"] = inv_cov

    return 0

def cleanup(config):
    pass
