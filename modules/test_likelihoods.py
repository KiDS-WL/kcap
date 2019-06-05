import numpy as np

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    likelihood_type = options.get_string(option_section, "likelihood", default="gaussian")
    likelihood_name = options.get_string(option_section, "likelihood_name", default=likelihood_type)
    n_dim = options.get_int(option_section, "n_dim", default=4)
    parameter_section = options.get_string(option_section, "parameter_section", default="params")
    return likelihood_type, n_dim, parameter_section, likelihood_name

def execute(block, config):
    likelihood_type, n_dim, parameter_section, likelihood_name = config
    theta = [block[parameter_section, f"param_{i}"] for i in range(n_dim)]
    
    mu = np.zeros(n_dim)
    inv_cov = np.eye(n_dim)
    d = theta-mu

    chi2 = np.einsum('i,ij,j', d, inv_cov, d)
    chi2 = float(chi2)
    like = -0.5*chi2

    block[names.data_vector, likelihood_name+"_CHI2"] = chi2
    block[names.likelihoods, likelihood_name+"_LIKE"] = like
    return 0

def clean(config):
    pass
