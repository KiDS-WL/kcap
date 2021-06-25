from __future__ import print_function
from cosmosis.datablock import names, option_section
import numpy as np
from scipy import special

def combonent(z, amp, z_mean, sigma):
    norm = np.sqrt(np.pi/2) * z_mean * sigma * special.erfc(-z_mean/(np.sqrt(2)*sigma)) + sigma**2 * np.exp(-z_mean**2/(2*sigma**2))
    return(np.exp(amp) * z * np.exp(-(z-z_mean)**2/(2*sigma**2)) / norm)
def comb_nz(z, amplitudes, means, sigma):
    return(np.sum([combonent(z, amplitudes[i], means[i], sigma) for i in range(len(amplitudes))], axis = 0))

def setup(options):
    amplitudes_file = options.get_string(option_section, "amplitudes_file")
    covariance_file = options.get_string(option_section, "covariance_file")
    z_min = options.get_double(option_section, "z_min")
    z_max = options.get_double(option_section, "z_max")
    sigma = options.get_double(option_section, "sigma")
    n_comp = options.get_int(option_section, "n_components")
    n_tomo = options.get_int(option_section, "n_tomo")
    nz_hist = options.get_int(option_section, "nz_hist")
    sum_amplitudes = options.get_bool(option_section, "sum_amplitudes")
    means = z_min + np.arange(n_comp)*sigma
    z = np.linspace(z_min, z_max, nz_hist)
    print("Loading amplitude data from {0}.".format(amplitudes_file))
    amplitudes = np.loadtxt(amplitudes_file)
    print("Loading covariance matrix of amplitude from {0}.".format(covariance_file))
    covariance = np.loadtxt(covariance_file)

    section = "NZ_SOURCE"
    data = {}
    if sum_amplitudes:
        print("Sum amplitudes to generate redshift distributions of {0} tomographic bins.".format(n_tomo))
        nz = np.zeros((n_tomo, len(z)))
        for i in range(n_tomo):
            nz[i,:] = comb_nz(z, amplitudes[i], means, sigma)
    else:
        print("Generate redshift distributions of {0} comb components.".format(n_comp))
        nz = np.zeros((n_comp, len(z)))
        for i in range(n_comp):
            nz[i,:] = combonent(z, 0, means[i], sigma)
    data[section] = (z, nz)
    return data, amplitudes, covariance, means, sigma


def execute(block, config):
    data_sets, amplitudes, covariance, means, sigma = config

    for name, data in list(data_sets.items()):
        z, nz = data
        nbin = len(nz)
        ns = len(z)
        block[name, "nbin"] = nbin
        block[name, "nz"] = ns
        block[name, "z"] = z
        for i, n in enumerate(nz):
            block[name, "bin_{0}".format(i + 1)] = n
    block["amplitudes", "amp"] = amplitudes
    block["amplitudes", "cov"] = covariance
    block["comb", "means"] = means
    block["comb", "sigma"] = sigma
    return 0
