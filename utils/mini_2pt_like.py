import collections

import scipy.interpolate
import numpy as np
pi = np.pi

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

def get_theory_point(x, y, mode="interpolate", interpolated_x=None, bin_edges=None, weighting=None, cut_bins=None):
    intp = scipy.interpolate.InterpolatedUnivariateSpline(np.log(x), y)
    
    if mode.lower() == "interpolate":
        result = intp(np.log(interpolated_x))
    elif mode.lower() == "integrate":
        result = np.zeros(bin_edges.size-1)
        for i in range(bin_edges.size-1):
            mask = (bin_edges[i] <= x) & (x < bin_edges[i+1])
            if isinstance(weighting, np.ndarray):
                w = weighting[mask]
                norm = np.trapz(w, x[mask])
            else:
                w = 1
                norm = 1
            result[i] = np.trapz(w*y[mask], x[mask])/norm
    else:
        raise ValueError(f"Mode {mode} not supported.")

    if cut_bins is not None:
        result = np.delete(result, cut_bins)

    return result




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

    angular_binning_mode = options.get_string(option_section, "angular_binning_mode", default="interpolate")
    if angular_binning_mode.lower() == "integrate":
        angular_bin_edges = options[option_section, "angular_bin_edges"]
        angular_bin_edges *= pi/180/60
    else:
        angular_bin_edges = None

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

    constant_c_offset = options.get_bool(option_section, "constant_c_offset", default=True)
    try:
        filename = options.get_string(option_section, "xi_pm_c_file")
        xi_p_c, xi_m_c = np.loadtxt(filename, unpack=True, usecols=[3,4])
        xi_p_c = np.delete(xi_p_c, cut_xi_plus_idx)
        xi_m_c = np.delete(xi_m_c, cut_xi_minus_idx)
    except errors.BlockNameNotFound:
        xi_p_c = None
        xi_m_c = None

    try:
        m_correction = options[option_section, "m_correction"]
        if len(m_correction) != n_bin:
            raise ValueError(f"Number of m-correction values supplied does not match number of tomographic bins: {len(m_correction)} vs {n_bin}.")
    except errors.BlockNameNotFound:
        m_correction = np.zeros(n_bin)

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

            data_vectors[i][j][0] *= 1/((1+m_correction[i])*(1+m_correction[j]))
            data_vectors[i][j][1] *= 1/((1+m_correction[i])*(1+m_correction[j]))

            data_vectors[i][j][0] = np.delete(data_vectors[i][j][0], cut_xi_plus_idx)
            data_vectors[i][j][1] = np.delete(data_vectors[i][j][1], cut_xi_minus_idx)
                
    like_name = options.get_string(option_section, "like_name")
    keep_theory_vector = options.get_bool(option_section, "keep_theory_vector", False)
    return inv_cov, data_vectors, theta_xi_plus, theta_xi_minus, \
           cut_xi_plus_idx, cut_xi_minus_idx, angular_binning_mode, angular_bin_edges, \
           order_cov, \
           constant_c_offset, xi_p_c, xi_m_c, like_name, keep_theory_vector

def execute(block, config):
    inv_cov, data_vectors, theta_xi_plus, theta_xi_minus, cut_xi_plus_idx, cut_xi_minus_idx, angular_binning_mode, angular_bin_edges, order_cov, constant_c_offset, xi_p_c, xi_m_c, like_name, keep_theory_vector = config
    
    n_bin = block["shear_xi_plus", "nbin_a"]

    if order_cov == "xi_pm-bin-theta":
        theory_xi_plus_vector = []
        theory_xi_minus_vector = []
        data_xi_plus_vector = []
        data_xi_minus_vector = []
        
        for i in range(n_bin):
            for j in range(i+1):
                data_xi_plus, data_xi_minus = data_vectors[i][j]

                theory_xi_plus_theta = block["shear_xi_plus", "theta"]
                theory_xi_plus = block["shear_xi_plus", f"bin_{i+1}_{j+1}"]
                theory_xi_minus_theta = block["shear_xi_minus", "theta"]
                theory_xi_minus = block["shear_xi_minus", f"bin_{i+1}_{j+1}"]

                theory_xi_plus = get_theory_point(theory_xi_plus_theta, theory_xi_plus,
                                                  mode=angular_binning_mode,
                                                  interpolated_x=theta_xi_plus,
                                                  bin_edges=angular_bin_edges,
                                                  weighting=theory_xi_plus_theta,
                                                  cut_bins=cut_xi_plus_idx)
                theory_xi_minus = get_theory_point(theory_xi_minus_theta, theory_xi_minus,
                                                   mode=angular_binning_mode,
                                                   interpolated_x=theta_xi_minus,
                                                   bin_edges=angular_bin_edges,
                                                   weighting=theory_xi_minus_theta,
                                                   cut_bins=cut_xi_minus_idx)

                if constant_c_offset:
                    theory_xi_plus += block["shear_c_bias", "delta_c"]
                    # theory_xi_minus += block["shear_c_bias", "delta_c"]
                if xi_p_c is not None:
                    theory_xi_plus += block["shear_c_bias", "A_c"]**2 * xi_p_c
                    theory_xi_minus += block["shear_c_bias", "A_c"]**2 * xi_m_c

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

                theory_xi_plus_theta = block["shear_xi_plus", "theta"]
                theory_xi_plus = block["shear_xi_plus", f"bin_{j+1}_{i+1}"]
                theory_xi_minus_theta = block["shear_xi_minus", "theta"]
                theory_xi_minus = block["shear_xi_minus", f"bin_{j+1}_{i+1}"]

                theory_xi_plus = get_theory_point(theory_xi_plus_theta, theory_xi_plus,
                                                  mode=angular_binning_mode,
                                                  interpolated_x=theta_xi_plus,
                                                  bin_edges=angular_bin_edges,
                                                  weighting=theory_xi_plus_theta,
                                                  cut_bins=cut_xi_plus_idx)
                theory_xi_minus = get_theory_point(theory_xi_minus_theta, theory_xi_minus,
                                                   mode=angular_binning_mode,
                                                   interpolated_x=theta_xi_minus,
                                                   bin_edges=angular_bin_edges,
                                                   weighting=theory_xi_minus_theta,
                                                   cut_bins=cut_xi_minus_idx)

                if constant_c_offset:
                    theory_xi_plus += block["shear_c_bias", "delta_c"]**2
                    # theory_xi_minus += block["shear_c_bias", "delta_c"]**2
                if xi_p_c is not None:
                    theory_xi_plus += block["shear_c_bias", "A_c"]**2 * xi_p_c
                    theory_xi_minus += block["shear_c_bias", "A_c"]**2 * xi_m_c

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
