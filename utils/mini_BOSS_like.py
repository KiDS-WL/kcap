import collections

import scipy.interpolate
import numpy as np
pi = np.pi

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors


def setup(options):
    data_file = options.get_string(option_section, "data_vector_file").split()
    
    data_vector = {}
    for i, f in enumerate(data_file):
        b = i + 1
        data_vector[b] = np.loadtxt(f)

        n_wedge = (data_vector[b].shape[1]-1)//2
        n_points_full = data_vector[b].shape[0]
        data_vector[b] = np.concatenate([data_vector[b][:,2*j+1] for j in range(n_wedge)])


    cov_file = options.get_string(option_section, "covariance_file").split()
    n_mocks = options.get_int(option_section, "n_mocks", default=2048)

    cov = {}
    for i, f in enumerate(cov_file):
        b = i + 1
        cov[b] = np.loadtxt(f)

    if set(data_vector.keys()) != set(cov.keys()):
        raise RuntimeError("Entries in data_file and covariance_file are not consistent.")
    try:
        points_range = options[option_section, "points_range"]
    except errors.BlockNameNotFound:
        points_range = [4, 32]
        
    cut_points_idx = list(range(points_range[0])) + list(range(points_range[1], n_points_full))

    if len(cut_points_idx) > 0:
        cut_points_idx = np.array(cut_points_idx, dtype=int)
        cut_points_idx = np.concatenate([cut_points_idx + n_points_full*i for i in range(n_wedge)])

        for b in data_vector.keys():
            data_vector[b] = np.delete(data_vector[b], cut_points_idx) 
            cov[b] = np.delete(cov[b], cut_points_idx, axis=0)
            cov[b] = np.delete(cov[b], cut_points_idx, axis=1)
    
    inv_cov = {b : np.linalg.inv(cov[b]) for b in cov.keys()}

    if n_mocks > 0:
        nu = n_mocks - 1
        for b, psi in inv_cov.items():
            p = psi.shape[0]
            inv_cov[b] = (nu-p-1)/nu * psi
    
    like_name = options.get_string(option_section, "like_name")
    keep_theory_vector = options.get_bool(option_section, "keep_theory_vector", False)
    return cov, inv_cov, data_vector, cut_points_idx, like_name, keep_theory_vector

def execute(block, config):
    cov, inv_cov, data_vector, cut_points_idx, like_name, keep_theory_vector = config
    
    # n_wedge = block["xi_wedges", "n_wedge"]
    # mu = block["xi_wedges", "vtheo_convolved"]
    chi2 = 0
    for b in data_vector.keys():
        mu = block[f"xi_wedges", f"bin_{b}"]
        mu = np.concatenate([m for m in mu])
        x = data_vector[b]
        d = x - mu

        chi2 += float(np.einsum('i,ij,j', d, inv_cov[b], d))

        if keep_theory_vector:
            x_sim = np.random.multivariate_normal(mu, cov[b])
            block[names.data_vector, like_name + f"_theory_bin_{b}"] = mu
            block[names.data_vector, like_name + f"_data_bin_{b}"] = x
            block[names.data_vector, like_name + f"_cov_bin_{b}"] = cov[b]
            block[names.data_vector, like_name + f"_inv_cov_bin_{b}"] = inv_cov[b]
            block[names.data_vector, like_name + f"_simulation_bin_{b}"] = x_sim
    
    like = -0.5*chi2

    block[names.data_vector, like_name+"_CHI2"] = chi2
    block[names.likelihoods, like_name+"_LIKE"] = like

    return 0

def cleanup(config):
    pass
