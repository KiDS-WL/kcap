import numpy as np

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

def setup(options):
    likelihood_type = options.get_string(option_section, "likelihood", default="gaussian")
    likelihood_name = options.get_string(option_section, "likelihood_name", default=likelihood_type)
    n_dim = options.get_int(option_section, "n_dim", default=4)
    parameter_section = options.get_string(option_section, "parameter_section", default="params")
    mu_section = options.get_string(option_section, "mu_section", default="mu")
    counter = [0,]
    return likelihood_type, n_dim, parameter_section, mu_section, likelihood_name, counter

def execute(block, config):
    likelihood_type, n_dim, parameter_section, mu_section, likelihood_name, counter = config
    theta = [block[parameter_section, f"param_{i}"] for i in range(n_dim)]
    
    mu = np.zeros(n_dim)
    # try:
    #     mu_slow = block[mu_section, f"mu_slow"]
    #     mu[:len(mu_slow)] = mu_slow
    # except errors.BlockError:
    #     pass

    inv_cov = np.eye(n_dim)
    d = theta-mu

    chi2 = np.einsum('i,ij,j', d, inv_cov, d)
    chi2 = float(chi2)
    like = -0.5*chi2

    block[names.data_vector, likelihood_name+"_CHI2"] = chi2
    block[names.likelihoods, likelihood_name+"_LIKE"] = like
    counter[0] += 1
    #print(counter[0])
    return 0

def clean(config):
    likelihood_type, n_dim, parameter_section, mu_section, likelihood_name, counter = config
    print(f"Ran likelihood {counter[0]} times.", flush=True)
